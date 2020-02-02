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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>

namespace OpenMS
{

  class OPENMS_DLLAPI SiriusMSFile
  {
public:

  // struct to store information about accsessions
  struct AccessionInfo
  {
    String sf_path;
    String sf_type;
    String sf_accession;
    String native_id_accession;
    String native_id_type;
  };

  // struct to store the compound information
  struct CompoundInfo
  {
    String cmp;
    double pmass;
    double rt;
    double fmz;
    String fid;
    String formula;
    int charge;
    String ionization;
    String des;
    String specref_format;
    String source_file;
    String source_format;
    std::vector<String> native_ids;
    std::vector<String> scan_indices;
    std::vector<String> specrefs;
  };

  /**
    @brief Internal structure used in @ref SiriusAdapter that is used
    for the conversion of a MzMlFile to an internal format.

    @ingroup ID

    Store .ms file.
    Comments (see CompoundInfo) are written to SIRIUS .ms file and additionally stores in CompoundInfo struct.
    If adduct information for a spectrum is missing, no adduct information is addded. 
    In this case, SIRIUS assumes default adducts for the respective spectrum.
    
    @return writes .ms file
    @return stores CompoundInfo
    
    @param spectra: Peakmap from input mzml.
    @param msfile: Writtes .ms file from sirius.
    @param feature_mapping: Adducts and features (index).
    @param feature_only: Only use features.
    @param isotope_pattern_iterations: At which depth to stop isotope_pattern extraction (if possible).
    @param v_cmpinfo: Vector of CompoundInfo.
    */

    // preprocessing e.g. feature information
    static void store(const PeakMap& spectra,
                      const OpenMS::String& msfile,
                      const FeatureMapping::FeatureToMs2Indices& feature_mapping,
                      const bool& feature_only,
                      const int& isotope_pattern_iterations,
                      const bool no_mt_info,
                      std::vector<SiriusMSFile::CompoundInfo>& v_cmpinfo);

  };

} // namespace OpenMS
