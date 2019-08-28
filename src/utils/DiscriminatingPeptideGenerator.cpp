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
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in_faa", "<file>", StringList(), "proteome files", true);
    setValidFormats_("in_faa", ListUtils::create<String>("fasta"));
    registerInputFile_("sample_sheet", "<file>", "", "sample sheet mapping filenames to samples and groups", true);
    setValidFormats_("in_back", ListUtils::create<String>("csv"));

    registerOutputFile_("out", "<file>", "", "group specific peptides");
    setValidFormats_("out", ListUtils::create<String>("csv"));

    registerOutputFile_("out_matrix", "<file>", "", "peptide abundance matrix");
    setValidFormats_("out_matrix", ListUtils::create<String>("csv"));

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


  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    const String markerfile_name = getStringOption_("in_marker");
    const String db_name = getStringOption_("in_back");
    const String outputfile_name_acc = getStringOption_("out");

    const Size min_size = static_cast<Size>(getIntOption_("min_length"));
    const Size max_size = static_cast<Size>(getIntOption_("max_length"));
    const Size missed_cleavages = static_cast<Size>(getIntOption_("missed_cleavages"));

    const unsigned isob = static_cast<Size>(getIntOption_("isob"));
    const unsigned target_tax_id = getIntOption_("target_tax_id");
    const String target_tax_level = getStringOption_("target_tax_level");

    //-------------------------------------------------------------
    // if given parse the nodes and names file
    //-------------------------------------------------------------
    TTaxMap nodes_map;
    TAccIdMap acc_id_map;
    const String nodesfile_name = getStringOption_("in_nodes");
    const String namesfile_name = getStringOption_("in_names");
    const String acc_id_filename = getStringOption_("in_acc_ids");

    bool use_ncbi = false;

    if (nodesfile_name != "" && namesfile_name != "" && acc_id_filename != "")
    {
        readNCBINodes(nodesfile_name, namesfile_name, target_tax_level, nodes_map);
        LOG_INFO << "DONE PARSING NODES" << std::endl;
        readIdMapping(acc_id_filename, acc_id_map);
        LOG_INFO << "DONE LOADING IDS" << std::endl;

        if (!nodes_map.count(target_tax_id))
        {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "target tax id not found",  std::to_string(target_tax_id));
        }
        use_ncbi = true;
    }

    //-------------------------------------------------------------
    // reading input marker peptides and build search index for Wu Manber
    //-------------------------------------------------------------
    std::vector<FASTAFile::FASTAEntry> marker_list;
    FASTAFile().load(markerfile_name, marker_list);

    seqan::String<seqan::CharString> needles;
#ifdef USE_DIGEST
    seqan::String<seqan::CharString> needles_dig;
#endif
    for (const auto & entry : marker_list)
    {
        std::string tmp = entry.sequence;
        if (isob > 0)
            std::replace(tmp.begin(), tmp.end(), 'I', 'L');
        if (isob > 1)
            std::replace(tmp.begin(), tmp.end(), 'Q', 'K');
#ifdef USE_DIGEST
        seqan::appendValue(needles_dig, seqan::CharString(std::string("!") + tmp + "#"));
#endif
        seqan::appendValue(needles, seqan::CharString(std::string(tmp)));
    }

    unsigned num_threads = 1;

#ifdef _OPENMP
#pragma omp parallel
{
    num_threads = omp_get_max_threads();
}
#endif

    typedef seqan::Pattern<seqan::String<seqan::CharString>, seqan::SetHorspool> TPattern;
    std::vector<std::shared_ptr<TPattern>> patterns_dig;

#ifdef USE_DIGEST
    typedef seqan::Pattern<seqan::String<seqan::CharString>, seqan::AhoCorasick> TPattern2;
    std::vector<std::shared_ptr<TPattern2>> patterns;
#endif

    for(size_t i = 0; i < num_threads; ++i)
    {
#ifdef USE_DIGEST
        patterns_dig.push_back(std::shared_ptr<TPattern>(new TPattern(needles_dig)));
#endif
        patterns.push_back(std::shared_ptr<TPattern2>(new TPattern2(needles)));
    }

    //-------------------------------------------------------------
    // iterate over database and search for tryptic peptides
    //-------------------------------------------------------------
    String enzyme = getStringOption_("enzyme");
    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme);
    digestor.setMissedCleavages(missed_cleavages);

    std::vector<std::list<String>> hits_digest(marker_list.size()), hits(marker_list.size());

    FASTAFile ff;
    ff.readStart(db_name);

    Size entry_nr = 0;
    const Size max_buffer_size = 1000000;
    std::vector<FASTAFile::FASTAEntry> buffer(max_buffer_size);
    while(!ff.atEnd())
    {
        //load fasta records into buffer for parallel processing
        buffer.clear();
        Size buff_size = 0;
        while(ff.readNext(buffer[buff_size]) && ++buff_size < max_buffer_size);

        const Size cbs = buff_size;
        //std::cerr << "TOP" << std::endl;
#ifdef _OPENMP
#pragma omp parallel
{
#endif
        std::vector<AASequence> current_digest;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for(Size i = 0; i < cbs; ++i)
        {
            const FASTAFile::FASTAEntry & fe = buffer[i];
#ifdef _OPENMP
            const unsigned thread_id = omp_get_thread_num();
#else
            const unsigned thread_id = 0;
#endif

#ifdef USE_DIGEST
            unsigned local_hit_cnt = 0, local_hit_cnt_dig = 0;
#endif
            seqan::CharString tmp(fe.sequence.c_str());
            if (isob > 0)
                std::replace(seqan::begin(tmp), seqan::end(tmp), 'I','L');
            if (isob > 1)
                std::replace(seqan::begin(tmp), seqan::end(tmp), 'Q','K');

            seqan::Finder<seqan::CharString> finder2(tmp);
//            seqan::Finder<seqan::CharString const> finder2(std::string(fe.sequence));
//            std::cerr << "Thread: " << thread_id << " : " <<  std::string(fe.sequence) << std::endl;
            //std::cerr << seqan::needle(*([thread_id])) << std::endl;
//            std::cerr << std::endl;
            while (seqan::find(finder2, *(patterns[thread_id])))
            {
//                std::cerr << "Hit2 " << seqan::beginPosition(finder2) << " : " << seqan::endPosition(finder2) << infix(finder2) << position(*(patterns[thread_id])) << std::endl;
//                std::cerr << "Hit2 " << fe.identifier << ": " << fe.sequence << " : " << infix(finder2) << std::endl;
                if(digestor.isValidProduct(fe.sequence, beginPosition(finder2), endPosition(finder2) - beginPosition(finder2), false, false, false))
                {
#ifdef _OPENMP
#pragma omp critical
#endif
//                    std::cerr << "Hit non" << fe.identifier << ": " << fe.sequence << " :: " << infix(finder2) << std::endl;
//                    std::cerr << "Hit non " << seqan::beginPosition(finder2) << " : " << seqan::endPosition(finder2) << " : " << infix(finder2) << std::endl;
//                    std::cerr << "Hit non " << fe.identifier << ": " << fe.sequence << " : " << infix(finder2) << std::endl;

                    hits[position(*(patterns[thread_id]))].push_back(fe.identifier);
#ifdef USE_DIGEST
                    ++local_hit_cnt;
#endif
                }
//                else
//                {
//                    std::cerr << "IVALID: " << fe.sequence << " " << beginPosition(finder2) << " " << endPosition(finder2) - beginPosition(finder2) << " " << infix(finder2) << std::endl;
//                    std::string tmp2 = fe.sequence.substr(beginPosition(finder2), endPosition(finder2) - beginPosition(finder2));
//                    std::cerr << tmp2 << std::endl;
//                    std::cerr << (std::find(current_digest.begin(), current_digest.end(), AASequence().fromString(tmp2)) != current_digest.end()) << std::endl;
//                }
            }
#ifdef USE_DIGEST
            if (enzyme == "none")
            {
                current_digest.push_back(AASequence::fromString(fe.sequence));
            }
            else
            {
                digestor.digest(AASequence::fromString(fe.sequence), current_digest, min_size, max_size);
            }

            for (const auto & seq : current_digest)
            {
                std::string tmp = seq.toString();
                if (isob > 0)
                    std::replace(tmp.begin(), tmp.end(), 'I', 'L');
                if (isob > 1)
                    std::replace(tmp.begin(), tmp.end(), 'Q', 'K');

                seqan::CharString haystack(seqan::CharString(std::string("!") + tmp + "#"));
                seqan::Finder<seqan::CharString> finder(haystack);

                //std::cout << "tid " <<  thread_id << " npats: " << patterns_dig.size() << std::endl;
                while (seqan::find(finder, *(patterns_dig[thread_id])))
                {
#ifdef _OPENMP
#pragma omp critical
#endif
//                    std::cerr << "Hit digest " << fe.identifier << ": " << fe.sequence << std::endl;
                    hits_digest[position(*(patterns_dig[thread_id]))].push_back(fe.identifier);
                    ++local_hit_cnt_dig;
                }
            }
/*            if (local_hit_cnt != local_hit_cnt_dig)
            {
                for (auto a : current_digest)
                    std::cout << a.toString() << std::endl;
                exit(1);
            }
*/
#endif
            ++entry_nr;
            if(!(entry_nr % 1000000))
                std::cout << entry_nr << '\n';
        }

#ifdef _OPENMP
}
#endif
        //std::cerr << "BOTTOM" << std::endl;
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------


    //original accessions
    String outputfile_name_tax(outputfile_name_acc);
    String outputfile_name_target(outputfile_name_acc);
    String outputfile_name_valid(outputfile_name_acc);
    outputfile_name_tax.substitute(".txt", "_tax.txt");
    outputfile_name_target.substitute(".txt", "_target.txt");
    outputfile_name_valid.substitute(".txt", "_valid.txt");

    TextFile out_acc, out_tax, out_target_tax, out_valid_hits;
    std::set<unsigned long> unique_tax, unique_target_tax;

    for(size_t i = 0; i < marker_list.size(); ++i)
    {
        out_acc << marker_list[i].sequence + ": " + ListUtils::concatenate(hits[i], ",");
        if(use_ncbi)
        {
            for (const auto & h : hits[i])
            {
                try
                {
                    unsigned long tax_id = getTaxIdFromFastaHeader(h, acc_id_map);
                    unique_tax.insert(tax_id);
                    const auto it = nodes_map.find(tax_id);
                    if (it != nodes_map.end())
                    {
                        if (nodes_map[tax_id].target_level_id == -1ul)
                        {
                            std::cerr << "INVALID TARGET LEVEL ID FOR " << tax_id << std::endl;
                            exit(1);
                        }

                        unique_target_tax.insert(nodes_map[tax_id].target_level_id);
                    }
                    else
                        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                            "No node entry found for tax id",  std::to_string(tax_id));
                } catch (Exception::InvalidValue & e)
                {
                    std::cerr << e.what() << "\n IGNORING FOR NOW!\n" << std::endl;
                }
            }
            std::vector<String> unique_tax_str, unique_target_tax_str;
            for(const auto & tid : unique_tax)
                unique_tax_str.push_back(nodes_map[tid].name);
            for(const auto & tid : unique_target_tax)
                unique_target_tax_str.push_back(nodes_map[tid].name);

            out_tax << marker_list[i].sequence + ": " + ListUtils::concatenate(unique_tax_str, ",");
            out_target_tax << marker_list[i].sequence + ": " + ListUtils::concatenate(unique_target_tax_str, ",");

            if (unique_target_tax.empty() || (unique_target_tax.size() == 1 && unique_target_tax.count(target_tax_id)))
                out_valid_hits << marker_list[i].sequence;

            unique_tax.clear();
            unique_target_tax.clear();
        }
    }
    out_acc.store(outputfile_name_acc);

    if(use_ncbi)
    {
        out_tax.store(outputfile_name_tax);
        out_target_tax.store(outputfile_name_target);
        out_valid_hits.store(outputfile_name_valid);
    }

#ifdef USE_DIGEST
    //-------------------------------------------------------------
    // Check for equality between normal and digest mode
    //-------------------------------------------------------------

    for(size_t i = 0; i < marker_list.size(); ++i)
    {
        bool eq = hits_digest[i].size() == hits[i].size();
        std::vector<String> tmp_hits( hits[i].begin(),  hits[i].end());
        std::vector<String> tmp_hits_digest( hits_digest[i].begin(),  hits_digest[i].end());
        std::sort(tmp_hits.begin(), tmp_hits.end());
        std::sort(tmp_hits_digest.begin(), tmp_hits_digest.end());

        eq = std::equal(tmp_hits_digest.begin(), tmp_hits_digest.end(), tmp_hits.begin());

        // std::list<String>::const_iterator it1(hits_digest[i].begin()), it2(hits[i].begin()), it1_end(hits_digest[i].end()), it2_end(hits[i].end());
        //
        // while (eq && it1 != it1_end && it2 != it2_end)
        // {
        //     eq = (*it1++ == *it2++);
        // }

        if (!eq)
            std::cerr << marker_list[i].sequence << ": " << hits_digest[i].size()
            << " " << hits[i].size() << " " << eq << std::endl;
    }
#endif

    return EXECUTION_OK;
  }


  struct SampleSheet{
      StringList filenames{};
      std::vector<unsigned> samples{};
      std::vector<unsigned> groups{};

      std::map<String, unsigned> filename_to_id{};
      std::map<String, unsigned> group_to_id{};
      std::map<String, unsigned> sample_to_id{};
  };

  void readSampleSheet(const String & filename, SampleSheet & samples)
  {
      //read the nodes - first col node id, sec col parent id, 3rd col rank
      std::ifstream is(filename);
      String line;
      StringList fields;
      Size rowCount = 0;

      TextFile::getLine(is, line); //skip header

      unsigned id = 0;
      while(TextFile::getLine(is, line))
      {

          line.split('\t', fields);

          if (fields.size() != 3)
              throw(Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "sample sheet line does not have 3 columns",  line));

          samples.filenames.push_back(fields[0]);
          samples.filename_to_id[fields[0]] = samples.filenames.size() - 1;

          auto it_s = samples.sample_to_id.find(fields[1]);
          if (it_s == samples.sample_to_id.end())
              std::tie(std::ignore, it_s) = samples.sample_to_id.insert(std::pair<String, unsigned>(fields[1], samples.sample_to_id.size()));

          samples.samples.push_back(it_s->second);

          auto it_g = samples.group_to_id.find(fields[2]);
          if (it_g == samples.group_to_id.end())
              std::tie(std::ignore, it_g) = samples.group_to_id.insert(std::pair<String, unsigned>(fields[2], samples.group_to_id.size()));

          samples.groups.push_back(it_g->second);
      }
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

  typedef std::map<String, std::vector<unsigned>> TPeptideMatrix;
  typedef std::map<String, std::vector<bool>> TPeptideBitMatrix;
  void generatePeptideMatrix(StringList faa, const SampleSheet & sheet, TPeptideMatrix & mat, const Options & opts)
  {
      ProteaseDigestion digestor;
      digestor.setEnzyme(opts.enzyme);
      digestor.setMissedCleavages(opts.missed_cleavages);

      const auto num_files = faa.size();
      const auto num_samples = sheet.sample_to_id.size();
      const auto num_groups = sheet.group_to_id.size();
      std::vector<unsigned> group_sizes(num_groups, 0);
      std::vector<unsigned> sample_sizes(num_samples, 0);

      TPeptideBitMatrix full_mat; //one column per file

      for (const String & faa : faa)
      {
          const auto file_id = sheet.filename_to_id.find(File::basename(faa))->second;
          const auto sample_id = sheet.samples[file_id];
          const auto group_id = sheet.groups[file_id];

          FASTAFile ff;
          ff.readStart(faa);

          FASTAFile::FASTAEntry fe;
          while(ff.readNext(fe))
          {
              std::vector<AASequence> current_digest;
              digestor.digest(AASequence::fromString(fe.sequence), current_digest, opts.min_size, opts.max_size);

              for (const auto & seq : current_digest)
              {
                  std::string tmp = seq.toString();
                  if (opts.isob > 0)
                      std::replace(tmp.begin(), tmp.end(), 'I', 'L');
                  if (opts.isob > 1)
                      std::replace(tmp.begin(), tmp.end(), 'Q', 'K');

//                  orig_pep_map[tmp].insert(seq.toString());
                  if (!full_mat.count(tmp))
                      full_mat[tmp] = std::vector<bool>(num_files, false);

                  full_mat[tmp][file_id] = true;
                  ++group_sizes[group_id];
                  ++sample_sizes[sample_id];
              }
          }
      }


      std::vector<unsigned> pep_cnt(num_samples, 0);
      std::vector<unsigned> group_cnt(num_groups, 0);
      for (const auto & row : full_mat)
      {
          //combine samples
          std::fill(pep_cnt.begin(), pep_cnt.end(), 0);
          for (Size i = 0; i < row.second.size(); ++i)
                pep_cnt[sheet.samples[i]] += row.second[i];

          for (Size i = 0; i < pep_cnt.size(); ++i)
              pep_cnt[i] = pep_cnt[i] >= opts.sample_t * sample_sizes[i] ? 1:0;

          //combine groups
          std::fill(group_cnt.begin(), group_cnt.end(), 0);
          for (Size i = 0; i < pep_cnt.size(); ++i)
              group_cnt[sheet.sample_to_group[i]] += pep_cnt[i];

          //is it group specific marker?
          bool above_t = false;
          unsigned pot_marker = -1;
          for (Size i = 0; i < group_cnt.size(); ++i)
          {
              if (group_cnt[i] >= opts.group_t * group_sizes[i])
              {
                  //already above threshold for another group?
                  if (above_t)
                  {
                      pot_marker = -1;
                      break;
                  }
                  pot_marker = i;
                  above_t = true;
              }
              else if (group_cnt[i] > 0) //below threshold but non-zero? KO!
              {
                  pot_marker = -1;
                  break;
              }
          }
          //store as marker
          if (pot_marker != -1)
              markers[pot_marker].push_back(row.first);
      }



      }

  }


int main(int argc, const char** argv)
{
  TOPPDiscriminatingPeptideGenerator tool;
  return tool.main(argc, argv);
}

/// @endcond
