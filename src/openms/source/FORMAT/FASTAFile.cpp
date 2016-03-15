// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Authors: Nico PFeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <fstream>

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

namespace OpenMS
{
  using namespace std;

  FASTAFile::FASTAFile()
  {

  }

  FASTAFile::~FASTAFile()
  {

  }

  void FASTAFile::load(const String& filename, vector<FASTAEntry>& data)
  {
    data.clear();

    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

    if (!File::readable(filename))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }


    seqan::SeqFileIn in(filename.c_str());
    //seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);

    std::string id_tmp, seq;
    String id, msg;
    String::size_type position = String::npos;
    Size size_read(0);

    while (!atEnd(in))
    {        
      try
      {
        readRecord(id_tmp, seq, in);

        if (data.empty()) msg = "The first entry could not be read!";
        else msg = "The last successful FASTA record was: '>" + data.back().identifier + "'. The record after failed.";

        FASTAEntry newEntry;
        newEntry.sequence = seq;
        newEntry.sequence.removeWhitespaces();

        // handle id
        id = id_tmp;
        id = id.trim();
        position = id.find_first_of(" \v\t");
        if (position == String::npos)
        {
          newEntry.identifier = id;
          newEntry.description = "";
        }
        else
        {
          newEntry.identifier = id.substr(0, position);
          newEntry.description = id.suffix(id.size() - position - 1);
        }
        id.clear();
        seq.clear();
        data.push_back(newEntry);
        size_read += newEntry.sequence.length();

      }
      catch (seqan::Exception const & e)
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", "Error while parsing FASTA file '" + filename + "'! " + msg +  " Please check the file!");
      }
    }

    if (size_read > 0 && data.empty())
      LOG_WARN << "No entries from FASTA file read. Does the file have MacOS "
               << "line endings? Convert to Unix or Windows line endings to"
               << " fix!" << std::endl;

    return;
  }

  void FASTAFile::store(const String& filename, const vector<FASTAEntry>& data) const
  {
    //ofstream outfile;
    //outfile.open(filename.c_str(), ofstream::out);
    seqan::SequenceOutputOptions fasta_opts;
    fasta_opts.lineLength = 80;

    seqan::SeqFileOut outfile;
    if(!seqan::open(outfile, filename.c_str()))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
    seqan::context(outfile).options.lineLength = 80;

    for (vector<FASTAEntry>::const_iterator it = data.begin(); it != data.end(); ++it)
    {      
      std::string header = it->identifier + " " + it->description;
      try
      {
        seqan::writeRecord(outfile, header, it->sequence.c_str());
      }
      catch(seqan::Exception & /*err*/)
      {
        throw Exception::IOException(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Error while writing FASTA file entry for id " + it->identifier + " to file " + filename);
      }
    }
  }
} // namespace OpenMS
