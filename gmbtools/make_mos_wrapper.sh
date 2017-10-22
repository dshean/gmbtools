#! /bin/bash

topdir=/nobackup/deshean/conus_combined/mos/conus_20171021_mos
cd $topdir

fn_list=""
fn_list+=" conus_mos_8m_summer_trans/conus_mos_8m_summer_trans.vrt"
fn_list+=" conus_mos_8m_latest_summer_trans/conus_mos_8m_latest_summer_trans.vrt"

gdalbuildvrt -r cubic conus_mos_8m_summmer.vrt $fn_list

fn_list=""
fn_list+=" conus_mos_8m/conus_mos_8m.vrt"
fn_list+=" conus_mos_8m_trans/conus_mos_8m_trans.vrt"
fn_list+=" conus_mos_8m_summer_trans/conus_mos_8m_summer_trans.vrt"
fn_list+=" conus_mos_8m_latest_summer_trans/conus_mos_8m_latest_summer_trans.vrt"

gdalbuildvrt -r cubic conus_mos_8m_all.vrt $fn_list




    
