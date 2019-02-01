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

#define USE_DIGEST

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MarkerPeptideChecker MarkerPeptideChecker

    @brief Checks for candidate marker peptides in background database.
<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MarkerPeptideChecker \f$ \longrightarrow \f$</td>
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
    @verbinclude UTILS_MarkerPeptideChecker.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MarkerPeptideCheckerr.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMarkerPeptideChecker :
  public TOPPBase
{
public:
  TOPPMarkerPeptideChecker() :
    TOPPBase("MarkerPeptideChecker", "Checks for candidate marker peptides in background database.", false)
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_marker", "<file>", "", "candidate peptides", true);
    setValidFormats_("in_marker", ListUtils::create<String>("fasta"));
    registerInputFile_("in_back", "<file>", "", "background DB", true);
    setValidFormats_("in_back", ListUtils::create<String>("fasta"));
    registerInputFile_("in_nodes", "<file>", "", "ncbi taxonomy nodes file", false);
    setValidFormats_("in_nodes", ListUtils::create<String>("csv"));
    registerInputFile_("in_names", "<file>", "", "ncbi taxonomy names file", false);
    setValidFormats_("in_names", ListUtils::create<String>("csv"));
    registerInputFile_("in_acc_ids", "<file>", "", "ncbi accessions to taxonomy ids file", false);
    setValidFormats_("in_acc_ids", ListUtils::create<String>("csv"));

    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("txt"));

    registerIntOption_("missed_cleavages", "<number>", 0, "The number of allowed missed cleavages", false);
    setMinInt_("missed_cleavages", 0);
    registerIntOption_("min_length", "<number>", 6, "Minimum length of peptide", false);
    registerIntOption_("max_length", "<number>", 40, "Maximum length of peptide", false);

    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("enzyme", "<string>", "Trypsin", "The type of digestion enzyme", false);
    setValidStrings_("enzyme", all_enzymes);    
  }

  struct NCBITax{
      std::string rank;
      std::string name;
      unsigned long parent_id;
      unsigned long target_level_id;
  };

  typedef std::map<unsigned long, NCBITax> TTaxMap;


  ExitCodes main_(int, const char**) override
  {

#ifdef _OPENMP
//    omp_set_dynamic(0);
//    omp_set_num_threads(16);

    std::cout << "#Threads: " << omp_get_num_threads() << " " << omp_get_max_threads() << std::endl;
#endif

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    const String markerfile_name = getStringOption_("in_marker");
    const String db_name = getStringOption_("in_back");
    const String outputfile_name_acc = getStringOption_("out");

    Size min_size = static_cast<Size>(getIntOption_("min_length"));
    Size max_size = static_cast<Size>(getIntOption_("max_length"));
    Size missed_cleavages = static_cast<Size>(getIntOption_("missed_cleavages"));

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
        readNCBINodes(nodesfile_name, namesfile_name, nodes_map);
        LOG_INFO << "DONE PARSING NODES" << std::endl;
        readIdMapping(acc_id_filename, acc_id_map);
        LOG_INFO << "DONE LOADING IDS" << std::endl;
        use_ncbi = true;
//        return ExitCodes::EXECUTION_OK;
    }

    //-------------------------------------------------------------
    // reading input marker peptides and build search index for Wu Manber
    //-------------------------------------------------------------

    std::vector<FASTAFile::FASTAEntry> marker_list;
    FASTAFile().load(markerfile_name, marker_list);

    seqan::String<seqan::CharString> needles, needles2;
    for (const auto & entry : marker_list)
    {
        seqan::appendValue(needles, seqan::CharString(std::string("!") + entry.sequence + "#"));
        seqan::appendValue(needles2, seqan::CharString(std::string(entry.sequence)));
    }

    unsigned num_threads = 1;

#ifdef _OPENMP
#pragma omp parallel
{
    num_threads = omp_get_max_threads();
}
#endif

    typedef seqan::Pattern<seqan::String<seqan::CharString>, seqan::SetHorspool> TPattern;
    typedef seqan::Pattern<seqan::String<seqan::CharString>, seqan::AhoCorasick> TPattern2;
    std::vector<std::shared_ptr<TPattern>> patterns;
    std::vector<std::shared_ptr<TPattern2>> patterns2;

    for(size_t i = 0; i < num_threads; ++i)
    {
        patterns.push_back(std::shared_ptr<TPattern>(new TPattern(needles)));
        patterns2.push_back(std::shared_ptr<TPattern2>(new TPattern2(needles2)));
    }
    //pattern(needles);


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

    unsigned entry = 0;
    unsigned max_buffer_size = 10000;
    std::vector<FASTAFile::FASTAEntry> buffer(max_buffer_size);
    while(!ff.atEnd())
    {
        //load fasta records into buffer for parallel processing
        buffer.clear();
        Size buff_size = 0;
        while(ff.readNext(buffer[buff_size]) && ++buff_size < max_buffer_size);

        const Size cbs = buff_size;
        std::cerr << "TOP" << std::endl;
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
//            std::cerr << i << std::endl;
            const FASTAFile::FASTAEntry & fe = buffer[i];

#ifdef _OPENMP
            const unsigned thread_id = omp_get_thread_num();
#else
            const unsigned thread_id = 0;
#endif
            seqan::CharString tmp = std::string(fe.sequence);
            seqan::Finder<seqan::CharString> finder2(tmp);
//            seqan::Finder<seqan::CharString const> finder2(std::string(fe.sequence));
//            std::cerr << "Thread: " << thread_id << " : " <<  std::string(fe.sequence) << std::endl;
            //std::cerr << seqan::needle(*(patterns2[thread_id])) << std::endl;
            while (seqan::find(finder2, *(patterns2[thread_id])))
            {
//                std::cerr << "Hit2 " << seqan::beginPosition(finder2) << " : " << seqan::endPosition(finder2) << infix(finder2) << position(*(patterns2[thread_id])) << std::endl;
//                std::cerr << "Hit2 " << fe.identifier << ": " << fe.sequence << " : " << infix(finder2) << std::endl;
                if(digestor.isValidProduct(fe.sequence, beginPosition(finder2), endPosition(finder2) - beginPosition(finder2), true, true, true))
                {
#ifdef _OPENMP
#pragma omp critical
#endif
                    hits[position(*(patterns2[thread_id]))].push_back(fe.identifier);
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
                seqan::CharString haystack(seqan::CharString(std::string("!") + seq.toString() + "#"));
                seqan::Finder<seqan::CharString> finder(haystack);

                //std::cout << "tid " <<  thread_id << " npats: " << patterns.size() << std::endl;
                while (seqan::find(finder, *(patterns[thread_id])))
                {
#ifdef _OPENMP
#pragma omp critical
#endif
//                    std::cerr << "Hit1 " << fe.identifier << ": " << fe.sequence << std::endl;
                    hits_digest[position(*(patterns[thread_id]))].push_back(fe.identifier);
                }
            }
#endif
//            ++entry;
//            if(!(entry % 1000))
//                std::cout << entry << '\n';
        }

#ifdef _OPENMP
}
#endif
        std::cerr << "BOTTOM" << std::endl;
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------


    //original accessions
    String outputfile_name_tax = outputfile_name_acc, outputfile_name_target = outputfile_name_acc;
    outputfile_name_tax.substitute(".txt", "_tax.txt");
    outputfile_name_target.substitute(".txt", "_target.txt");

    TextFile out_acc, out_tax, out_target_tax;
    std::set<unsigned long> unique_tax, unique_target_tax;

    for(size_t i = 0; i < marker_list.size(); ++i)
    {
        out_acc << marker_list[i].sequence + ": " + ListUtils::concatenate(hits[i], ",");
        if(use_ncbi)
        {
            for (const auto & h : hits_digest[i])
            {
                try
                {
                    unsigned long tax_id = getTaxIdFromFastaHeader(h, acc_id_map);
                    unique_tax.insert(tax_id);
                    unique_target_tax.insert(nodes_map[tax_id].target_level_id);
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
            unique_tax.clear();
            unique_target_tax.clear();
        }
    }
    out_acc.store(outputfile_name_acc);

    if(use_ncbi)
    {
        out_tax.store(outputfile_name_tax);
        out_target_tax.store(outputfile_name_target);
    }

#ifdef USE_DIGEST
    //-------------------------------------------------------------
    // DEBUG OUTPUT
    //-------------------------------------------------------------

    for(size_t i = 0; i < marker_list.size(); ++i)
    {
        bool eq = hits_digest[i].size() == hits[i].size();

        std::list<String>::const_iterator it1(hits_digest[i].begin()), it2(hits[i].begin()), it1_end(hits_digest[i].end()), it2_end(hits[i].end());

        while (eq && it1 != it1_end && it2 != it2_end)
        {
            eq = (*it1++ == *it2++);
        }

        if (!eq)
            std::cout << marker_list[i].sequence << ": " << hits_digest[i].size() << " " << hits[i].size() << " " << eq << std::endl;
    }
#endif

    return EXECUTION_OK;
  }


  typedef std::unordered_map<std::string, unsigned long> TAccIdMap;

  inline unsigned long getTaxIdFromFastaHeader(const String & fasta_header, const TAccIdMap & id_map)
  {
      const auto finder = id_map.find(extractNcbiAccession(fasta_header));
      if (finder != id_map.end())
          return finder->second;
      else
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              String("No tax id found for accession [") + extractNcbiAccession(fasta_header) + "]", fasta_header);
  }

  inline String extractNcbiAccession(const String & fasta_header)
  {
      StringList tmp;
      fasta_header.split('|', tmp);
      if (tmp.size() > 1)
          return tmp[1];
      else
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Invalid fasta header", fasta_header);
  }

  void readIdMapping(const String & filename, TAccIdMap & map)
  {
      //read the nodes - first col node id, sec col parent id, 3rd col rank
      std::ifstream is(filename);
      String line;
      StringList fields;
      Size rowCount = 0;
      while(TextFile::getLine(is, line))
      {
          line.split('\t', fields);
          map[fields[0]] = stoul(fields[1]);
#ifdef DEBUG_MARKER
          if (!(++rowCount % 1000000))
              std::cerr << rowCount << std::endl;
#endif
      }
  }

  void readNCBINodes(const String & filename_nodes,
                     const String & filename_names,
                     TTaxMap & nodes_map)
  {
      //read the nodes - first col node id, sec col parent id, 3rd col rank
      CsvFile csv_nodes(filename_nodes, '|');

      //read the nodes
      for (size_t i = 0; i < csv_nodes.rowCount(); ++i)
      {
          StringList row;
          csv_nodes.getRow(i, row);
          nodes_map[stoul(row[0])] = {row[2], "", stoul(row[1]), -1ul};
      }

      //set ids for target level (species)
      for(auto & node : nodes_map)
      {
          unsigned long next_id = node.first;
          unsigned long target_id = -1ul;
          std::vector<unsigned long> tmp_ids;

          while(target_id == -1ul)
          {
              if (nodes_map[next_id].target_level_id != -1ul)
              {
                  target_id = nodes_map[next_id].target_level_id;
                  break;
              }
              tmp_ids.push_back(next_id);
              if (nodes_map[next_id].rank == "species")
                  target_id = next_id;
              else if (nodes_map[next_id].parent_id != -1ul && nodes_map[next_id].parent_id != next_id)
                  next_id = nodes_map[next_id].parent_id;
              else
                  break;
          }
          if (target_id != -1ul)
          {
              for(const auto id : tmp_ids)
              {
                  nodes_map[id].target_level_id = target_id;
              }
          }
      }

      //get the names for ids
      CsvFile csv_names(filename_names, '|');
      for (size_t i = 0; i < csv_names.rowCount(); ++i)
      {
          StringList row;
          csv_names.getRow(i, row);
          if (nodes_map[stoul(row[0])].name.empty() || row[3] == "scientific name")
              nodes_map[stoul(row[0])].name = row[1];
      }

#ifdef DEBUG_MARKER
      //DEBUG output
      for(const auto & node : nodes_map)
      {
        std::cout << node.first << ": " << node.second.name << " " << node.second.rank;
        if (node.second.target_level_id != -1ul)
            std::cout << " -> (" << node.second.target_level_id << ") " << nodes_map[node.second.target_level_id].name;
        std::cout << std::endl;
      }
#endif

  }

};


int main(int argc, const char** argv)
{
  TOPPMarkerPeptideChecker tool;
  return tool.main(argc, argv);
}

/// @endcond
