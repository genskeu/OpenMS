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
// $Maintainer: Ulrich Genske $
// $Authors: Ulrich Genske $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCRTPeakWidth.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

//Destructor
QCRTPeakWidth::~QCRTPeakWidth()
{

}

//Constructor
QCRTPeakWidth::QCRTPeakWidth(std::vector<OpenMS::FeatureMap> files):
  maps(files)
  {

  }

//Main method to write mztab peptide section data needed for MBR alignment plot (PTXQC)
int QCRTPeakWidth::RTPeakWidth(MzTab& mztab) const
{

  //vector<FeatureMap> maps;
  MzTabPSMSectionRows rows;
  MzTabPSMSectionRows mztabRowsPSM = mztab.getPSMSectionRows();
  vector<MzTabString> unique_ids_;

  //maps = feat_map_;

  for (Size m = 0; m < maps.size(); m++)
  {

    String rfile;

    if (maps[m].metaValueExists("spectra_data"))
    {
      StringList rfiles = maps[m].getMetaValue("spectra_data");
      rfile = rfiles[0];
    }

    for (vector<Feature>::const_iterator f_it = maps[m].begin(); f_it!=maps[m].end();f_it++)
    {

      MzTabPSMSectionRow row;

      //set column retention time
      MzTabDouble retention_time;
      retention_time.set(  f_it->getRT() );
      vector<MzTabDouble> retention_times;
      retention_times.push_back(retention_time);

      MzTabDoubleList retention_time_mztab;
      retention_time_mztab.set(retention_times);
      row.retention_time = retention_time_mztab;

      //set column mass to charge
      MzTabDouble exp_mz;
      exp_mz.set(f_it->getMZ());
      row.exp_mass_to_charge = exp_mz;


    	//Set optional columns: RT length, unique_id, source file and intensity
      vector<MzTabOptionalColumnEntry> v_opt;

      double retention_length;
      retention_length = f_it->getMetaValue("FWHM");
      MzTabOptionalColumnEntry RTlength = make_pair("opt_retention_length",MzTabString(retention_length));
      v_opt.push_back (RTlength);

      UInt64 unique_id = f_it->getUniqueId();
      MzTabOptionalColumnEntry u_id = make_pair("opt_unique_id",MzTabString(unique_id));
      v_opt.push_back (u_id);

      MzTabOptionalColumnEntry sraw = make_pair("opt_raw_source_file",MzTabString(rfile));
      v_opt.push_back (sraw);

      double intensity;
      intensity = f_it->getIntensity();
      MzTabOptionalColumnEntry opt_intensity = make_pair("opt_intensity",MzTabString(intensity));
      v_opt.push_back (opt_intensity);

      row.opt_ = v_opt;


      rows.push_back(row);


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
        newMZTabrows.push_back(*it_new_row);
        continue;
      }

      bool inMzTab = FALSE;

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

        //if id matches combine rows
        if (UID_mzTab.toCellString().compare(UID_new.toCellString())==0)
        {
          inMzTab = TRUE;
          //add new data from the new row to the exisiting one only its empty
          if(it_mzTab_row->sequence.isNull()){it_mzTab_row->sequence = it_new_row->sequence;}
          if(it_mzTab_row->retention_time.isNull()){it_mzTab_row->retention_time = it_new_row->retention_time;}
          if(it_mzTab_row->charge.isNull()){it_mzTab_row->charge = it_new_row->charge;}
          if(it_mzTab_row->spectra_ref.isNull()){it_mzTab_row->spectra_ref = it_new_row->spectra_ref;}
          if(it_mzTab_row->exp_mass_to_charge.isNull()){it_mzTab_row->exp_mass_to_charge = it_new_row->exp_mass_to_charge;}

          //add all optional columns (does not check for dublicates)
          for (vector<MzTabOptionalColumnEntry>::const_iterator o_it = it_new_row->opt_.begin(); o_it != it_new_row->opt_.end(); ++o_it)
          {
            it_mzTab_row->opt_.push_back(*o_it);
          }
        }
      }
      //add new rows with no match
      if(inMzTab == FALSE)
      {
        //no match
        newMZTabrows.push_back(*it_new_row);
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
