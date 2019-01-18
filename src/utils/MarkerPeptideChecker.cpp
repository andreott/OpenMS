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
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

#include<seqan/find.h>

#include <map>

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
    registerInputFile_("in_acc_tid_map", "<file>", "", "map accessions to taxonomic id", false);
    setValidFormats_("in_acc_tid_map", ListUtils::create<String>("csv"));
    registerInputFile_("in_tid_tname_map", "<file>", "", "map taxonomic ids to taxonomic names", false);
    setValidFormats_("in_tid_tname_map", ListUtils::create<String>("csv"));

    registerOutputFile_("out", "<file>", "", "Output file (peptides)");
    setValidFormats_("out", ListUtils::create<String>("fasta"));

    registerIntOption_("missed_cleavages", "<number>", 1, "The number of allowed missed cleavages", false);
    setMinInt_("missed_cleavages", 0);
    registerIntOption_("min_length", "<number>", 6, "Minimum length of peptide", false);
    registerIntOption_("max_length", "<number>", 40, "Maximum length of peptide", false);

    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("enzyme", "<string>", "Trypsin", "The type of digestion enzyme", false);
    setValidStrings_("enzyme", all_enzymes);    
  }


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
    String markerfile_name = getStringOption_("in_marker");
    String db_name = getStringOption_("in_back");
    String outputfile_name = getStringOption_("out");

    Size min_size = getIntOption_("min_length");
    Size max_size = getIntOption_("max_length");
    Size missed_cleavages = getIntOption_("missed_cleavages");

    //-------------------------------------------------------------
    // reading input marker peptides and build search index for Wu Manber
    //-------------------------------------------------------------

    std::vector<FASTAFile::FASTAEntry> marker_list;
    FASTAFile().load(markerfile_name, marker_list);

    seqan::String<seqan::CharString> needles;
    for (const auto & entry : marker_list)
    {
        seqan::appendValue(needles, seqan::CharString(std::string("!") + entry.sequence + "#"));
    }

    unsigned num_threads = 1;

#ifdef _OPENMP
#pragma omp parallel
{
    num_threads = omp_get_max_threads();
}
#endif

    typedef seqan::Pattern<seqan::String<seqan::CharString>, seqan::SetHorspool> TPattern;
    //TPattern p(needles);
    std::vector<std::shared_ptr<TPattern>> patterns;

    for(size_t i = 0; i < num_threads; ++i)
    {
//        std::shared_ptr<TPattern> pp(new TPattern(needles));
        patterns.push_back(std::shared_ptr<TPattern>(new TPattern(needles)));
    }
    //pattern(needles);


    //-------------------------------------------------------------
    // iterate over database and search for tryptic peptides
    //-------------------------------------------------------------
    String enzyme = getStringOption_("enzyme");
    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme);
    digestor.setMissedCleavages(missed_cleavages);

    std::vector<std::list<String>> hits(marker_list.size());

    FASTAFile ff;
    ff.readStart(db_name);

    unsigned entry = 0;
    unsigned max_buffer_size = 10000;
    std::vector<FASTAFile::FASTAEntry> buffer(max_buffer_size);
    while(!ff.atEnd())
    {
        //load fasta records into buffer for parallel processing
        buffer.clear();
        unsigned buff_size = 0;
        while(ff.readNext(buffer[buff_size]) && ++buff_size < max_buffer_size);

        const int cbs = buff_size;
        std::cerr << "TOP" << std::endl;
#ifdef _OPENMP
#pragma omp parallel
{
#endif
        std::vector<AASequence> current_digest;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for(int i = 0; i < cbs; ++i)
        {
            const FASTAFile::FASTAEntry & fe = buffer[i];

            if (enzyme == "none")
            {
                current_digest.push_back(AASequence::fromString(fe.sequence));
            }
            else
            {
                digestor.digest(AASequence::fromString(fe.sequence), current_digest, min_size, max_size);
            }
#ifdef _OPENMP
            const unsigned thread_id = omp_get_thread_num();
#else
            const unsigned thread_id = 1;
#endif
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
                    hits[position(*(patterns[thread_id]))].push_back(fe.identifier);
                }
            }

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


    for(size_t i = 0; i < marker_list.size(); ++i)
    {
        std::cout << marker_list[i].sequence << ": ";
        for (const auto & h : hits[i])
        {
            std::cout << h << ",";
        }
        std::cout << '\n';
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPMarkerPeptideChecker tool;
  return tool.main(argc, argv);
}

/// @endcond
