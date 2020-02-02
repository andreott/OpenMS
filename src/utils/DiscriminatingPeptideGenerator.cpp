// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Sandro Andreotti $
// $Authors: Sandro Andreotti $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

#include<seqan/find.h>

#include <map>
#include <unordered_map>
#include <algorithm>

#include <QDir>

#define USE_DIGEST

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_DiscriminatingPeptideGenerator DiscriminatingPeptideGenerator

    @brief Identifies group specific peptides from a list of amino acid fasta files.
<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ DiscriminatingPeptideGenerator \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none (FASTA input) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none (FASTA input) none (FASTA output) </td>
        </tr>
    </table>
</CENTER>

    This application is used to check a list of peptides against a background database (probably huge)
    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_DiscriminatingPeptideGenerator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_DiscriminatingPeptideGeneratorr.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDiscriminatingPeptideGenerator :
  public TOPPBase
{
public:
  TOPPDiscriminatingPeptideGenerator() :
    TOPPBase("DiscriminatingPeptideGenerator", "Identifies group specific peptides from a list of amino acid fasta files.", false)
  {

  }

protected:
  typedef std::map<String, std::vector<unsigned>> TPeptideMatrix;
  //typedef std::map<String, std::vector<bool>> TPeptideBitMatrix;
  typedef std::unordered_map<std::string, std::vector<bool>> TPeptideBitMatrix;
  typedef std::unordered_map<std::string, std::vector<Size>> TPeptideProteinTie;
  typedef std::vector<std::vector<std::pair<String, StringList>>> TMarkerList; //one list of peptides per group and one list of assiciated protein accessions
  typedef std::unordered_map<std::string, std::set<std::string>> TPepMap; //stores mapping from isobaric modified peptide to original pep sequences

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in_faa", "<file>", StringList(), "proteome files", true);
    setValidFormats_("in_faa", ListUtils::create<String>("fasta"));
    registerInputFile_("sample_sheet", "<file>", "", "sample sheet mapping filenames to samples and groups", true);
    setValidFormats_("sample_sheet", ListUtils::create<String>("csv"));

    registerOutputFile_("out", "<file>", "", "group specific peptides");
    setValidFormats_("out", ListUtils::create<String>("csv"));

    registerStringOption_("out_faa_dir", "<file>", "", "if set for each input file a file with the same name containing markers will be stored to this directory.", false, false);

    registerOutputFile_("out_matrix", "<file>", "", "peptide abundance matrix");
    setValidFormats_("out_matrix", ListUtils::create<String>("csv"));
    
    registerOutputFile_("out_mapping", "<file>", "", "original to marker peptide mapping", false);
    setValidFormats_("out_mapping", ListUtils::create<String>("csv"));

    registerIntOption_("missed_cleavages", "<number>", 0, "The number of allowed missed cleavages", false);
    setMinInt_("missed_cleavages", 0);
    registerIntOption_("min_length", "<number>", 6, "Minimum length of peptide", false);
    registerIntOption_("max_length", "<number>", 40, "Maximum length of peptide", false);

    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("enzyme", "<string>", "Trypsin", "The type of digestion enzyme", false);
    setValidStrings_("enzyme", all_enzymes);

    registerIntOption_("isob", "<number>", 0, "account for isobarics: 0 - none, 1 I/L, 2 I/L + K/Q", false);
    setMinInt_("isob", 0);

    registerDoubleOption_("sample_threshold", "<double>", 0.6, "fraction of lists for the same sample that must contain a peptide", false);
    setMinFloat_("sample_threshold", 0.);
    setMaxFloat_("sample_threshold", 1.);

    registerDoubleOption_("group_threshold", "<double>", 1.0, "fraction of lists for the same group that must contain a peptide", false);
    setMinFloat_("group_threshold", 0.);
    setMaxFloat_("group_threshold", 1.);
  }

  struct Options
  {
      Size min_size{};
      Size max_size{};
      Size isob{};
      Size missed_cleavages{};
      String enzyme{};
      String output_dir{};
      double sample_t{};
      double group_t{};
  };

  ExitCodes main_(int, const char**) override
  {

    Options opt;
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    const StringList faa_files =  getStringList_("in_faa");
    const String samsheet_file = getStringOption_("sample_sheet");
    const String outfile = getStringOption_("out");
    const String outfile_mapping = getStringOption_("out_mapping");
    const String out_dir = getStringOption_("out_faa_dir");

    if (!out_dir.empty())
    {
        QDir dir(out_dir.toQString());

        if (!dir.exists())
        {
            LOG_ERROR << "Specified path does not exist\n";
            return ExitCodes::CANNOT_WRITE_OUTPUT_FILE;
        }
        opt.output_dir = dir.absolutePath();
    }

    opt.min_size = static_cast<Size>(getIntOption_("min_length"));
    opt.max_size = static_cast<Size>(getIntOption_("max_length"));
    opt.enzyme = getStringOption_("enzyme");
    opt.missed_cleavages = static_cast<Size>(getIntOption_("missed_cleavages"));
    opt.isob = static_cast<Size>(getIntOption_("isob"));
    opt.sample_t = getDoubleOption_("sample_threshold");
    opt.group_t = getDoubleOption_("group_threshold");

    //-------------------------------------------------------------
    // reading sample sheet and validate against input faa file list
    //-------------------------------------------------------------

    SampleSheet sheet;
    if(!readSampleSheet(samsheet_file, faa_files, sheet))
        return ExitCodes::INPUT_FILE_CORRUPT;
    std::cerr << "done reading sample sheet\n";

    //-------------------------------------------------------------
    // generate marker peptides
    //-------------------------------------------------------------

    TPeptideMatrix mat;
    TMarkerList marker;
    TPepMap orig_pep_map;
    generatePeptideMatrix(faa_files, sheet, mat, opt, marker, orig_pep_map);
    std::cerr << "done generating markers\n";

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    TextFile of;
    FASTAFile of_fasta;
    for (Size i = 0; i < sheet.group_names.size(); ++i)
    {
        if (!marker[i].empty())
        {
            StringList seqs(marker[i].size());
            std::transform(marker[i].begin(), marker[i].end(), seqs.begin(), [](decltype (marker[i][0]) & p){return p.first;});
            of << sheet.group_names[i] + '\t' + ListUtils::concatenate(seqs,"\t");
        }
        else
            of << sheet.group_names[i];

        if (!opt.output_dir.empty())
        {
            std::vector<FASTAFile::FASTAEntry> fa_entries;
            for (const auto & m : marker[i])
                fa_entries.push_back(FASTAFile::FASTAEntry(m.first,ListUtils::concatenate(m.second, ' '), m.first));
            of_fasta.store(QDir(opt.output_dir.toQString()).filePath(sheet.group_names[i].toQString().replace(' ','_') + ".fa"), fa_entries);
        }
    }
    of.store(outfile);
    
    if (!outfile_mapping.empty())
    {
        TextFile ofm;
        for (const auto & entry : orig_pep_map)
        {            
            ofm << entry.first + "," + ListUtils::concatenate(entry.second, ',');
        }
        ofm.store(outfile_mapping);
    }


    return EXECUTION_OK;
  }


  struct SampleSheet{
      StringList filenames{}; //filenames
      std::vector<unsigned> samples{}; //file -> sample
      std::vector<unsigned> groups{}; //sample -> group

      std::map<String, unsigned> filename_to_id{};
      StringList sample_names;
      StringList group_names;
  };

  bool readSampleSheet(const String & filename, const StringList filenames_faa, SampleSheet & samples)
  {
      //read the nodes - first col node id, sec col parent id, 3rd col rank
      std::ifstream is(filename);
      String line;
      StringList fields;

      StringList basenames(filenames_faa.size());
      std::transform(filenames_faa.begin(), filenames_faa.end(), basenames.begin(), File::basename);
      std::sort(basenames.begin(), basenames.end());

      std::map<String, unsigned> group_to_id{}, sample_to_id{};
      TextFile::getLine(is, line); //skip header

      while(TextFile::getLine(is, line))
      {          
          line.split('\t', fields);

          if (fields.size() != 3)
              throw(Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "sample sheet line does not have 3 columns",  line));

          if (!ListUtils::contains(basenames, fields[0]))
          {
              LOG_WARN << "Not all filenames in sample sheet were given as input! Missing Filename: " << fields[0] << '\n';
              continue;
          }

          samples.filenames.push_back(fields[0]);
          samples.filename_to_id[fields[0]] = samples.filenames.size() - 1;

          //check if group was already seen. if not new id
          auto it_g = group_to_id.find(fields[2]);
          if (it_g == group_to_id.end())
          {
              samples.group_names.push_back(fields[2]);
              std::tie(it_g, std::ignore) = group_to_id.insert(std::pair<String, unsigned>(fields[2], group_to_id.size()));
          }

          auto it_s = sample_to_id.find(fields[1]);
          if (it_s == sample_to_id.end()) //new sample -> add sample -> group connection
          {
              samples.sample_names.push_back(fields[1]);
              std::tie(it_s, std::ignore) = sample_to_id.insert(std::pair<String, unsigned>(fields[1], sample_to_id.size()));
              samples.groups.push_back(it_g->second);
          }
          else
          {
              //sanity check that this sample has not previously been assigned to other group
              if (samples.groups[it_s->second] != it_g->second)
              {
                  throw(Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "ambiguous sample -> group asignment!",  ""));
              }
          }
          samples.samples.push_back(it_s->second);
      }

      //OH OH - some filenames have no entry in sample sheet! Find out and report their names and return error
      if (samples.filenames.size() < basenames.size())
      {
          StringList diff_a, fn_sorted(samples.filenames);
          std::sort(fn_sorted.begin(), fn_sorted.end());
          std::back_insert_iterator<StringList> it(diff_a);
          std::set_difference(basenames.begin(), basenames.end(), samples.filenames.begin(), samples.filenames.end(), it);

          LOG_ERROR << "Not all filenames in sample sheet!\nMissingFilenames: " << ListUtils::concatenate(diff_a, "\t");
          return false;
      }
      return true;
  }

  bool validateInput(StringList faa, const SampleSheet & sheet)
  {
      StringList basenames(faa.size());
      std::transform(faa.begin(), faa.end(), basenames.begin(), File::basename);
      std::sort(basenames.begin(), basenames.end());

      StringList sheet_fn = sheet.filenames;
      std::sort(sheet_fn.begin(), sheet_fn.end());

      StringList diff_a;
      std::back_insert_iterator<StringList> it(diff_a);
      std::set_difference(basenames.begin(), basenames.end(), sheet_fn.begin(), sheet_fn.end(), it);

      if(!diff_a.empty())
      {
          LOG_ERROR << "Not all filenames in sample sheet!\nMissingFilenames: " << ListUtils::concatenate(diff_a, "\t");
          return false;
      }

      std::set_difference(sheet_fn.begin(), sheet_fn.end(), basenames.begin(), basenames.end(), it);
      if(!diff_a.empty())
      {
          LOG_WARN << "Not all filenames in sample sheet were given as input!\nMissingFilenames: " << ListUtils::concatenate(diff_a, "\t");
      }

      return true;
  }


  void generatePeptideMatrix(StringList faa_list,
                             const SampleSheet & sheet,
                             TPeptideMatrix & mat,
                             const Options & opts,
                             TMarkerList & markers,
                             TPepMap & orig_pep_map)
  {
      ProteaseDigestion digestor;
      digestor.setEnzyme(opts.enzyme);
      digestor.setMissedCleavages(opts.missed_cleavages);

      const auto num_files = faa_list.size();
      const auto num_samples = sheet.sample_names.size();
      const auto num_groups = sheet.group_names.size();
      std::vector<unsigned> group_sizes(num_groups, 0);
      std::vector<unsigned> sample_sizes(num_samples, 0);

      std::vector<String> accessions; //stores protein accessions
      TPeptideProteinTie acc_trace; //match peptides to protein accessions

      TPepMap orig_pep_map_full;
      TPeptideBitMatrix full_mat; //one column per file
      std::vector<bool> sample_check(num_samples, false);

      for (const String & faa : faa_list)
      {
          const auto file_id = sheet.filename_to_id.find(File::basename(faa))->second;
          const auto sample_id = sheet.samples[file_id];
          ++sample_sizes[sample_id];

          std::cerr << sample_id << " " << sheet.samples[file_id] << std::endl;

          if(!sample_check[sample_id])
          {
              sample_check[sample_id] = true;
              ++group_sizes[sheet.groups[sample_id]];
          }

          FASTAFile ff;
          ff.readStart(faa);

          FASTAFile::FASTAEntry fe;
          while(ff.readNext(fe))
          {

              accessions.push_back(fe.identifier);
              std::vector<AASequence> current_digest;

              //remove marker for stop codon as the created 'X' cause trouble and would have to be deleted from peptides
              if(fe.sequence.back() == '*')
                  fe.sequence.pop_back();

              digestor.digest(AASequence::fromString(fe.sequence, false), current_digest, opts.min_size, opts.max_size);

              for (const auto & seq : current_digest)
              {
                  std::string tmp = seq.toString();
                  if (opts.isob > 0)
                      std::replace(tmp.begin(), tmp.end(), 'I', 'L');
                  if (opts.isob > 1)
                      std::replace(tmp.begin(), tmp.end(), 'Q', 'K');

                  orig_pep_map_full[tmp].insert(seq.toString());
                  if (!full_mat.count(tmp))                  
                      full_mat[tmp] = std::vector<bool>(num_files, false);

                  full_mat[tmp][file_id] = true;
                  acc_trace[tmp].push_back(accessions.size()-1);
              }
          }
      }

      std::vector<unsigned> sample_cnt(num_samples, 0), sample_cnt_bin(num_samples, 0);
      std::vector<unsigned> group_cnt_cons(num_groups, 0), group_cnt(num_groups, 0);

      markers.clear();
      markers.resize(num_groups);

      const Size undef = static_cast<Size>(-1);
      for (const auto & row : full_mat)
      {
          //combine samples
          std::fill(sample_cnt.begin(), sample_cnt.end(), 0);
          std::fill(sample_cnt_bin.begin(), sample_cnt_bin.end(), 0);
          for (Size i = 0; i < row.second.size(); ++i)
                sample_cnt[sheet.samples[i]] += row.second[i];

          for (Size i = 0; i < sample_cnt.size(); ++i)
              sample_cnt_bin[i] = sample_cnt[i] >= opts.sample_t * sample_sizes[i] ? 1:0;

          //combine groups
          std::fill(group_cnt.begin(), group_cnt.end(), 0);
          std::fill(group_cnt_cons.begin(), group_cnt_cons.end(), 0);
          for (Size i = 0; i < sample_cnt.size(); ++i)
          {
              group_cnt[sheet.groups[i]] += sample_cnt[i];
              group_cnt_cons[sheet.groups[i]] += sample_cnt_bin[i];
          }

          //is it group specific marker?
          bool above_t = false;
          Size pot_marker = undef;
          for (Size i = 0; i < group_cnt.size(); ++i)
          {
              if (group_cnt_cons[i] >= opts.group_t * group_sizes[i])
              {
                  //already above threshold for another group?
                  if (above_t)
                  {
                      pot_marker = undef;
                      break;
                  }
                  pot_marker = i;
                  above_t = true;
              }
              else if (group_cnt[i] > 0) //below threshold but non-zero for some replicate? KO!
              {
                  pot_marker = undef;
                  break;
              }
          }
          //store as marker
          if (pot_marker != static_cast<Size>(-1))
          {
              StringList tmp_acc;
              for (const auto a : acc_trace[row.first])
                  tmp_acc.push_back(accessions[a]);

              markers[pot_marker].push_back(std::make_pair(row.first, std::move(tmp_acc)));
              orig_pep_map.insert({row.first, orig_pep_map_full[row.first]});
          }
      }
  }
};


int main(int argc, const char** argv)
{
  TOPPDiscriminatingPeptideGenerator tool;
  return tool.main(argc, argv);
}

/// @endcond