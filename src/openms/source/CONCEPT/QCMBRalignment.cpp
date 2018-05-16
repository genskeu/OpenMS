// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Maria Trofimova $
// $Authors: Maria Trofimova $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCMBRalignment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

//Destructor
QCMBRalignment::~QCMBRalignment()
{

}

//Constructor
QCMBRalignment::QCMBRalignment(std::vector<OpenMS::FeatureMap> files):
  maps(files)
  {

  }

//Main method to write mztab peptide section data needed for MBR alignment plot (PTXQC)
int QCMBRalignment::MBRAlignment(MzTab& mztab) const
{

  //vector<FeatureMap> maps;
  MzTabPSMSectionRows rows;
  MzTabPSMSectionRows mztabRowsPSM = mztab.getPSMSectionRows();
  vector<MzTabString> unique_ids_;
  int pepIDCount = 0;

  //maps = feat_map_;

  for (Size m = 0; m < maps.size(); m++)
  {

    String rfile;
    //Keep index of current map for spectra reference
    Size run = m;

    if (maps[m].metaValueExists("spectra_data"))
    {
      StringList rfiles = maps[m].getMetaValue("spectra_data");
      rfile = rfiles[0];
    }

    for (vector<Feature>::const_iterator f_it = maps[m].begin(); f_it!=maps[m].end();f_it++)
    {

      vector<PeptideIdentification> pep_id = f_it->getPeptideIdentifications();

      if (pep_id.empty())
      {
        //If we want to write empty lines peptide_data_ for features with retention times
        //
      }

      else
      {
        //Iterate over peptide hits of a feature and write data into mztab
        //In case of 2 and more hits take the first one
 	    for (vector<PeptideIdentification>::iterator p_it = pep_id.begin(); p_it!=pep_id.end(); p_it++)
 	    {
 	      pepIDCount++;
 	      MzTabPSMSectionRow row;
 	      MzTabString PepSeq;
          MzTabDouble correctRT;
          MzTabDouble oriRT;

 		  //Set sequence
          vector<PeptideHit> hits = p_it->getHits();
          PeptideHit hit = hits[0];
          AASequence seq = hit.getSequence();
          PepSeq.set(seq.toString());
          row.sequence = PepSeq;

          //Set corrected RTs
          float corrRT = f_it->getRT();
          correctRT.set(corrRT);
          vector<MzTabDouble> cRTs;
          cRTs.push_back (correctRT);
          MzTabDoubleList listC;
          listC.set(cRTs);
          row.retention_time = listC;

    	  float orRT;
    	  if (p_it->metaValueExists("original_RT")) {orRT = p_it->getMetaValue("original_RT");}
 	      else {throw "Retention time was not written in the feature.";}

 	      //Set spectrum reference
    	  String spectrum_ref = p_it->getMetaValue("spectrum_reference");
    	  const String& const_spec = spectrum_ref;
    	  MzTabSpectraRef spec_ref;
    	  Size index = run+1;
    	  spec_ref.setMSFile(index);
    	  spec_ref.setSpecRef(spectrum_ref);
    	  spec_ref.setSpecRefFile(const_spec);
    	  row.spectra_ref = spec_ref;

    	  //Set optional columns: original RT and source file
    	  String ori = to_string(orRT);
          MzTabOptionalColumnEntry oRT = make_pair("opt_original_retention_time",MzTabString(ori));
          vector<MzTabOptionalColumnEntry> v;
          v.push_back (oRT);

          UInt64 id = f_it->getUniqueId();
          MzTabOptionalColumnEntry u_id = make_pair("opt_unique_id",MzTabString(id));
          v.push_back (u_id);
          unique_ids_.push_back (MzTabString(id));

          MzTabOptionalColumnEntry sraw = make_pair("opt_raw_source_file",MzTabString(rfile));
          v.push_back (sraw);
          row.opt_ = v;

          rows.push_back(row);

         }
 	  }

 	}

  }


  //Merge new lines and existing lines. Based on unique ids (UniqueIdInterface)
  //If PSM section was not written before: append rows

  if (mztabRowsPSM.empty())
  {
    mztab.setPSMSectionRows(rows);
  }

  //Else: append rows
  //If unique ids are equal: append columns (relevant for metric)
  //If not: insert unique id in dictionary and add new line to mztab
  else
  {
    // vector saving the rows which are completly new (no id match)
    vector<MzTabPSMSectionRow> newMZTabrows;

    for (vector<MzTabPSMSectionRow>::iterator it_mzTab_row = mztabRowsPSM.begin(); it_mzTab_row != mztabRowsPSM.end(); ++it_mzTab_row)
    {
      //check unique id of mzTab row
      vector<MzTabOptionalColumnEntry> opt_row_mzTab = it_mzTab_row->opt_;
      MzTabString UID_mzTab;

      for (vector<MzTabOptionalColumnEntry>::const_iterator o_it = opt_row_mzTab.begin(); o_it != opt_row_mzTab.end();++o_it)
      {
        if(o_it->first=="opt_unique_id")
          {
            UID_mzTab = o_it->second;
          }
      }
      //if there is no id continue
      if(UID_mzTab.isNull())
      {
        continue;
      }

      bool inMzTab = FALSE;

      for (vector<MzTabPSMSectionRow>::const_iterator it_new_row = rows.begin(); it_new_row != rows.end(); ++it_new_row)
      {
        //check unique id of new rows
        vector<MzTabOptionalColumnEntry> opt_row_new = it_new_row->opt_;
        MzTabString UID_new;


        for (vector<MzTabOptionalColumnEntry>::const_iterator o_it = opt_row_new.begin(); o_it != opt_row_new.end();++o_it)
        {
          if(o_it->first=="opt_unique_id")
            {
              UID_new = o_it->second;
            }
        }
        //if there is no id in the new row add row and continue
        if(UID_new.isNull())
        {
          //only add rows with no id once
          if(it_mzTab_row == mztabRowsPSM.begin()) newMZTabrows.push_back(*it_new_row);
          continue;
        }



        //if id matches combine rows
        if (UID_mzTab.toCellString().compare(UID_new.toCellString())==0)
        {
          inMzTab = TRUE;
          //add new data from the new row to the exisiting one
          if(it_mzTab_row->sequence.isNull()){it_mzTab_row->sequence = it_new_row->sequence;}
          if(it_mzTab_row->retention_time.isNull()){it_mzTab_row->retention_time = it_new_row->retention_time;}
          if(it_mzTab_row->charge.isNull()){it_mzTab_row->charge = it_new_row->charge;}
          if(it_mzTab_row->spectra_ref.isNull()){it_mzTab_row->spectra_ref = it_new_row->spectra_ref;}
          //add all optional columns (does not check for dublicates)
          for (vector<MzTabOptionalColumnEntry>::const_iterator o_it = it_new_row->opt_.begin(); o_it != it_new_row->opt_.end(); ++o_it)
          {
            it_mzTab_row->opt_.push_back(*o_it);
          }
        }

        //add new rows with no match
        if(inMzTab == FALSE && it_new_row == rows.end() )
        {
          //no match
          newMZTabrows.push_back(*it_new_row);
        }
      }
    }
    //combine the updated old rows and the new ones
    vector<MzTabPSMSectionRow> mztabRowsPSMnew = mztabRowsPSM;
    mztabRowsPSMnew.insert(mztabRowsPSMnew.end(), newMZTabrows.begin() , newMZTabrows.end() );
    //setPSMSection
    mztab.setPSMSectionRows(mztabRowsPSMnew);
  }

  return rows.size()!=0 ? 1:0;
}
