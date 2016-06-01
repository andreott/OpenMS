// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXDELTAMASSESGENERATOR_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXDELTAMASSESGENERATOR_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief generates complete list of all possible mass shifts due to isotopic labelling
   * 
   * Isotopic labelling results in the shift of peptide masses.
   * 
   * For example in a Lys8/Arg10 SILAC labelled sample, some peptides (the ones with one
   * Arg in their sequence) will show a relative mass shift between light and heavy
   * partners of 10 Da. This class constructs the complete list of all possible mass
   * shifts that arise from isotopic labelling.
   */
  class OPENMS_DLLAPI MultiplexDeltaMassesGenerator
  {
    public:

    /**
     * @brief constructor
     * 
     * @param labels    isotopic labels, e.g. '[][Lys8,Arg10]' for SILAC labelling. For a detailed description and
     * further examples see the algorithm::labels parameter in the FeatureFinderMultiplex und MultiplexResolver tools.
     * @param missed_cleavages    maximum number of missed cleavages due to incomplete digestion
     * @param label_mass_shift    name of labels and their corresponding mass shifts
     * For example due to knock-outs in one of the samples.
     */
    MultiplexDeltaMassesGenerator(String labels, int missed_cleavages, std::map<String,double> label_mass_shift);
        
    /**
     * @brief generate all mass shifts that can occur due to the absence of one or multiple peptides
     * (e.g. for a triplet experiment generate the doublets and singlets that might be present)
     */
    void generateKnockoutDeltaMasses();

    /**
     * @brief write the list of labels for each of the sample
     */
    void printSamplesLabelsList() const;
    
    /**
     * @brief write the list of all mass patterns
     */
    void printDeltaMassesList() const;
    
    /**
     * @brief returns the list of mass shift patterns
     */
    std::vector<MultiplexDeltaMasses> getDeltaMassesList();
    
    /**
     * @brief returns the list of mass shift patterns
     */
    const std::vector<MultiplexDeltaMasses>& getDeltaMassesList() const;
    
    /**
     * @brief returns the list of samples with their corresponding labels
     */
    std::vector<std::vector<String> > getSamplesLabelsList();
    
    /**
     * @brief returns the list of samples with their corresponding labels
     */
    const std::vector<std::vector<String> >& getSamplesLabelsList() const;
    
    /**
     * @brief returns the short label string
     * 
     * @param label    long label, e.g. "Label:13C(6)15N(4)"
     */
    String getLabelShort(String label);
    
    /**
     * @brief returns the long label string
     * 
     * @param label    short label, e.g. "Arg10"
     */
    String getLabelLong(String label);
    
    /**
     * @brief extract the label set from the sequence
     *
     * @param sequence    amino acid sequence
     */
    MultiplexDeltaMasses::LabelSet extractLabelSet(AASequence sequence);
    
    private:
   
    /**
     * @brief isotopic labels
     */
    String labels_;
    
    /**
     * @brief flat list of all occuring isotopic labels
     */
    std::vector<String> labels_list_;
    
    /**
     * @brief list of samples with their corresponding labels
     */
    std::vector<std::vector<String> > samples_labels_;
    
    /**
     * @brief maximum number of missed cleavages
     */
    int missed_cleavages_;

    /**
     * @brief list of all possible mass shift patterns
     */
    std::vector<MultiplexDeltaMasses> delta_masses_list_;
      
    /**
     * @brief mapping from single label to delta mass
     * e.g. "Arg10" -> 10.0082686
     */
    std::map<String,double> label_delta_mass_;
    
    /**
     * @brief mapping from a short label (as in the user params) to a long label (as in PSI-MS name)
     * e.g. "Arg10" -> "Label:13C(6)15N(4)"
     */
    std::map<String,String> label_short_long_;
    
    /**
     * @brief mapping from a long label (as in PSI-MS name) to a short label (as in the user params)
     * e.g. "Label:13C(6)15N(4)" -> "Arg10"
     */
    std::map<String,String> label_long_short_;
    
 };
  
}

#endif /* MULTIPLEXDELTAMASSESGENERATOR_H */
