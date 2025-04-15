;
; read original CEDS/CMIP6 SSP emissions files
; regrid to 1 deg and convert from kg/m2/s to molecules/cm2/s
; sum all anthro sectors to one variable
;  (requires gnu compiler, not intel, on cheyenne/geyser)
;

begin
dir_inp = "/glade/p/acom/MUSICA/emissions/cams/CAMS-GLOB-ANTv4.2/orig/"
version= "v4.2"
;dir_out = "/glade/u/home/tilmes/acom-climate/tilmes/emis/CAMS-GLOB-ANTv4.2/2000_2020/new/"
dir_out = "/glade/scratch/$USER/CAMS-GLOB-ANTv4.2/2000_2020/"
resolution = "0.9x1.25"
 PI  = get_pi("float" )

print(dir_out)
print(dir_inp)
;
; constants
;
  avog = 6.022e23
  earth_rad  = 6.371e8 ; in cm
  earth_area = 4. * acos(-1.) * earth_rad * earth_rad
; define creation date
;
      creation_date = systemfunc ("date +%Y%m%d")
;

;
; template for lat-lon CESM grid
;
  in = addfile("/glade/p/cesmdata/cseg/inputdata/atm/cam/inic/fv/template_" + resolution + "_L32.nc","r")
  lat_out = in->lat
  lon_out = in->lon
  nlat_out = dimsizes(lat_out)
  nlon_out = dimsizes(lon_out)
  delete(in)
;
; create weights from lat (to be used for global average)
;
  gw = NormCosWgtGlobe(lat_out)
;
; make sure output dir exists
;
  system("mkdir -p " + dir_out)
;
; loop over species
;
   do nspc=0,30
;
    if ( nspc.eq.0 ) then
      spc_name = "BC"
      spc_name1 = "bc"
      spc_name2 = "bc_a4"
      mw = 12.
      sf = 1.0
    end if
    if ( nspc.eq.1 ) then
      spc_name = "CO"
      spc_name1 = "co"
      spc_name2 = "CO"
      mw = 28.
      sf = 1.0
    end if
    if ( nspc.eq.2 ) then
      spc_name = "NH3"
      spc_name1 = "nh3"
      spc_name2 = "NH3"
      mw = 17.
      sf = 1.0
    end if
    if ( nspc.eq.3 ) then
      spc_name = "NOx"
      spc_name1 = "nox"
      spc_name2 = "NOx"
      mw = 30.
      ;sf = 1.0
      ; different factor required because of mistake in original emissions
      sf = 46./30. 
    end if
    if ( nspc.eq.4 ) then
      spc_name = "OC"
      spc_name1 = "oc"
      spc_name2 = "pom_a4"
       sf = 1.4
      mw = 12.
    end if
    if ( nspc.eq.5 ) then
      spc_name = "SO2"
      spc_name1 = "so2"
      spc_name2 = "SO2"
      mw = 64.
      sf = 1.0
    end if
    if ( nspc.eq.6 ) then
     spc_name = "C2H6"
     spc_name1 = (/"ethane"/)
     mw = 30.
     sf = 1.0
    end if
    if ( nspc.eq.7 ) then
     spc_name = "C3H8"
     spc_name1 = (/"propane"/)
     mw = 44.
     sf = 1.0
   end if
   if ( nspc.eq.8 ) then
     spc_name = "BIGALK"
      spc_name1 = (/"butanes","pentanes","hexanes","esters","ethers"/)
      mw = 72.
      mw_spec_set = (/58., 72., 86., 184., 81./) 
      sf = 1.0
   end if
    if ( nspc.eq.9 ) then
     spc_name = "C2H4"
      spc_name1= (/"ethene"/)
      mw = 28.
      sf = 1.0
   end if
    if ( nspc.eq.10 ) then
     spc_name = "C3H6"
      spc_name1 = (/"propene"/)
      mw = 42.
      sf = 1.0
   end if
    if ( nspc.eq.11 ) then
     spc_name = "C2H2"
      spc_name1 = (/"ethyne"/)
      mw = 26.
      sf = 1.0
   end if
    if ( nspc.eq.12) then
     spc_name = "BIGENE"
      ;spc_name1 = (/"olef"/)
      spc_name1 = (/"other-alkenes-and-alkynes"/)
      mw = 56.
      sf = 1.0
   end if
   if ( nspc.eq.13 ) then
     spc_name = "BENZENE"
      spc_name1= (/"benzene"/)
      mw = 78.
      sf = 1.0
   end if
    if ( nspc.eq.14 ) then
     spc_name = "TOLUENE"
      spc_name1 = (/"toluene"/)
      mw = 92.
      sf = 1.0
   end if
    if ( nspc.eq.15 ) then
     spc_name = "XYLENES"
      spc_name1 = (/"xylene","trimethylbenzene","other-aromatics"/)
      mw = 106.
      mw_spec_set = (/106., 120., 126./) 
      sf = 1.0
   end if
    if ( nspc.eq.16 ) then
     spc_name = "CH2O"
      spc_name1 = (/"formaldehyde"/)
      mw = 30.
      sf = 1.0
   end if
    if ( nspc.eq.17 ) then
     spc_name = "CH3CHO"
      spc_name1 = (/"other-aldehydes"/)
      mw = 44.
      sf = 1.
   end if
  if ( nspc.eq.18 ) then
     spc_name = "CH3OH"
      spc_name1 = (/"alcohols"/)
      sf = 0.15
      mw = 32.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.19 ) then
     spc_name = "C2H5OH"
      spc_name1 = (/"alcohols"/)
      sf = 0.85
      mw = 46.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.20 ) then
     spc_name = "CH3COCH3"
      spc_name1 = (/"ketones"/)
      sf = 0.2
      mw = 58.
   end if
    if ( nspc.eq.21 ) then
     spc_name = "MEK"
      spc_name1 = (/"ketones"/)
      sf = 0.8
      mw = 72.
   end if
    if ( nspc.eq.22 ) then
     spc_name = "HCOOH"
      spc_name1 = (/"acids"/)
      sf = 0.5
      mw = 46.
   end if
  if ( nspc.eq.23 ) then
     spc_name = "CH3COOH"
      spc_name1 = (/"acids"/)
      sf = 0.5
      mw = 60.
   end if
    if ( nspc.eq.24 ) then
     spc_name = "ISOP"
      spc_name1 = (/"isoprene"/)
      sf = 1.
      mw = 68.
   end if
    if ( nspc.eq.25 ) then
     spc_name = "MTERP"
      spc_name1  =(/"monoterpenes"/)
      sf = 1.
      mw = 136.
   end if
    if ( nspc.eq.26 ) then
     spc_name = "HCN"
      spc_name1 = "co"
      sf = 0.003
      mw = 27.
    end if
    if ( nspc.eq.27 ) then
     spc_name = "CH3CN"
     spc_name1 = "co"
      sf = 0.002
      mw = 41.
    end if
    if ( nspc.eq.28) then
     spc_name = "IVOC"
     spc_name1 = "IVOC"
     sf = 1.
     mw = 184.
    end if
    if ( nspc.eq.29) then
     spc_name = "SVOC"
     spc_name1 = "SVOC"
     sf = 1.
     mw = 310.
    end if
    if ( nspc.eq.30) then
     spc_name = "SOAGx1.5"
     spc_name1 = "SOAGx1.5"
     sf = 1.5
     mw = 12.
    end if

    nspec1 = dimsizes(spc_name1)
    str1 = str_join(spc_name1,"+")
    if sf .lt. 1 then
       str2 = sprintf("%4.2e", sf)
    else
      str2= " "
    end if
    hist_spec = "CAMS species: "+str2+str1
    print(hist_spec)

if nspc .le. 27 then

  do ced  = 0, nspec1-1
;
; get anthro files
    files = systemfunc("ls " + dir_inp + "CAMS-GLOB-ANT_2000-2020_0.1x0.1_anthro_"+spc_name1(ced) + ".nc")
    nfils = dimsizes(files)
    print(files)
; should be just 1 file per species
    nfil=0
;
      file_inp = files(nfil)
      in = addfile(file_inp,"r")
   if ced.eq.0 then
      lat_inp   = in->lat
      lon_inp   = in->lon
      time      = in->time
      ntim = dimsizes(time)
      time0 = ispan(0,ntim-1,1)
      time0 = (/ 1315200, 1315944, 1316640, 1317384, 1318104, 1318848, 1319568,  \
    1320312, 1321056, 1321776, 1322520, 1323240, 1323984, 1324728, 1325400,  \
    1326144, 1326864, 1327608, 1328328, 1329072, 1329816, 1330536, 1331280,  \
    1332000, 1332744, 1333488, 1334160, 1334904, 1335624, 1336368, 1337088,  \
    1337832, 1338576, 1339296, 1340040, 1340760, 1341504, 1342248, 1342920,  \
    1343664, 1344384, 1345128, 1345848, 1346592, 1347336, 1348056, 1348800,  \
    1349520, 1350264, 1351008, 1351704, 1352448, 1353168, 1353912, 1354632,  \
    1355376, 1356120, 1356840, 1357584, 1358304, 1359048, 1359792, 1360464,  \
    1361208, 1361928, 1362672, 1363392, 1364136, 1364880, 1365600, 1366344,  \
    1367064, 1367808, 1368552, 1369224, 1369968, 1370688, 1371432, 1372152,  \
    1372896, 1373640, 1374360, 1375104, 1375824, 1376568, 1377312, 1377984,  \
    1378728, 1379448, 1380192, 1380912, 1381656, 1382400, 1383120, 1383864,  \
    1384584, 1385328, 1386072, 1386768, 1387512, 1388232, 1388976, 1389696,  \
    1390440, 1391184, 1391904, 1392648, 1393368, 1394112, 1394856, 1395528,  \
    1396272, 1396992, 1397736, 1398456, 1399200, 1399944, 1400664, 1401408,  \
    1402128, 1402872, 1403616, 1404288, 1405032, 1405752, 1406496, 1407216,  \
    1407960, 1408704, 1409424, 1410168, 1410888, 1411632, 1412376, 1413048,  \
    1413792, 1414512, 1415256, 1415976, 1416720, 1417464, 1418184, 1418928,  \
    1419648, 1420392, 1421136, 1421832, 1422576, 1423296, 1424040, 1424760,  \ 
    1425504, 1426248, 1426968, 1427712, 1428432, 1429176, 1429920, 1430592,  \
    1431336, 1432056, 1432800, 1433520, 1434264, 1435008, 1435728, 1436472,  \
    1437192, 1437936, 1438680, 1439352, 1440096, 1440816, 1441560, 1442280,  \
    1443024, 1443768, 1444488, 1445232, 1445952, 1446696, 1447440, 1448112,  \
    1448856, 1449576, 1450320, 1451040, 1451784, 1452528, 1453248, 1453992,  \
    1454712, 1455456, 1456200, 1456896, 1457640, 1458360, 1459104, 1459824,  \
    1460568, 1461312, 1462032, 1462776, 1463496, 1464240, 1464984, 1465656,  \
    1466400, 1467120, 1467864, 1468584, 1469328, 1470072, 1470792, 1471536,  \
    1472256, 1473000, 1473744, 1474416, 1475160, 1475880, 1476624, 1477344,  \
    1478088, 1478832, 1479552, 1480296, 1481016, 1481760, 1482504, 1483176,  \
    1483920, 1484640, 1485384, 1486104, 1486848, 1487592, 1488312, 1489056,  \
    1489776, 1490520, 1491264, 1491960, 1492704, 1493424, 1494168, 1494888,  \
    1495632, 1496376, 1497096, 1497840, 1498560 /) 
      
      delete(time)
      time=time0
      time!0="time"
      time@units="hours since 1850-01-01 00:00:00"
      time@calendar="gregorian"
     printVarSummary(time)
;
; define date array
;
     
      utc_date = cd_calendar(time, 0)
      date = doubletointeger(time)
      date@units = "YYYYMMDD"
      date@long_name = "date"
      date!0 = "time"
      date&time = time
      date = tointeger(utc_date(:,0)) * 10000 + tointeger(utc_date(:,1)) * 100 + tointeger(utc_date(:,2))
      print(date)
  end if
;

;
      print("Processing " + file_inp)
;
; find emissions variable name
;
      vnames = getfilevarnames(in)
      nnames = dimsizes(vnames)
print(vnames)

; get sum of all 
   var1 = in->ind
   mw_spec = stringtofloat(var1@molecular_weight)
   print(mw_spec)
   if mw_spec.ne.mw then
     print("mw_spec different from mw")
   end if
   mw_s = mw
   if nspec1 .gt. 1 then
     mw_s = mw_spec_set(ced)
   end if

;
;
; convert from kg/m2/s to molec/cm2/s
;
;                            kg to g         m2 to cm2
;                               |                |
;     emiss_inp = emiss_inp * 1.e3 * avog/mw_s* 1.e-4

;if min(emiss_inp).lt.0  then
;     print("negative data")
;     wk1D      = ndtooned(emiss_inp)
;     ibad = ind(wk1D.lt.0)
;     wk1D(ibad)=0.
;     emiss_inp = onedtond(wk1D, dimsizes(emiss_inp))
; else
;     print("no negative data")
; end if


;
; regrid to CESM lat-lon grid
;
;     emiss_out = area_conserve_remap_Wrap(lon_inp,lat_inp,emiss_out1,lon_out,lat_out,False)
;     delete(emiss_out1)
;     delete(emiss_out@units)
;     emiss_out@units = "molecules/cm2/s"
;     emiss_out@molecular_weight = mw_spec
;     emiss_out@molecular_weight_units = "g/mole"


;
; compute global average (in molec/cm2/s)
; annual average in Tg/year
;
;     emiss_out_gl = wgt_areaave(emiss_out,gw,1.0,1)
;     dims = dimsizes(emiss_out_gl)
;     nyr = dims(0)/12
;     emiss_out_gl_ann = emiss_out_gl(0:nyr-1)
;     do n=0,nyr-1
;       n1 = n * 12
;       n2 = n1 + 11
;       emiss_out_gl_ann(n) = dim_avg_Wrap(emiss_out_gl(n1:n2))
;     end do
;     emiss_out_gl_ann = emiss_out_gl_ann * earth_area * 86400. * 365 * mw_spec * 1.e-12/avog
;    print(emiss_out_gl_ann)

;     delete(emiss_out_gl)
;     delete(emiss_out_gl_ann)



    print("read ene")
      sect_energy = in->ene 
    print("read ind")
      sect_industry = in->ind 
      sect_industry2 = in->fef 
      sect_solvents= in->slv
      sect_transport = in->tro 
      sect_transport2 = in->tnr 
      sect_resident = in->res 
      sect_agriculture= in->awb
      sect_agriculture2= in->ags
      sect_agriculture3= in->agl
      sect_waste = in->swd 
    print("read shp")
      sect_shipping = in->shp 
    
    printVarSummary(sect_energy)
    ; fillv = sect_energy@_FillValue

    ; sect_agriculture=where(sect_agriculture.eq.fillv,0.,sect_agriculture)
    ; sect_energy=where(sect_energy.eq.fillv,0.,sect_energy)
    ; sect_industry=where(sect_industry.eq.fillv,0.,sect_industry)
    ; sect_industry2=where(sect_industry2.eq.fillv,0.,sect_industry2)
    ; sect_solvents=where(sect_solvents.eq.fillv,0.,sect_solvents)
    ; sect_transport=where(sect_transport.eq.fillv,0.,sect_transport)
    ; sect_transport2=where(sect_transport2.eq.fillv,0.,sect_transport2)
    ; sect_resident=where(sect_resident.eq.fillv,0.,sect_resident)
    ; sect_agriculture2=where(sect_agriculture2.eq.fillv,0.,sect_agriculture2)
    ; sect_agriculture3=where(sect_agriculture3.eq.fillv,0.,sect_agriculture3)
    ; sect_waste=where(sect_waste.eq.fillv,0.,sect_waste)
    ; sect_shipping=where(sect_shipping.eq.fillv,0.,sect_shipping)
    ; delete(fillv)

   if ced.eq.0 then
      sect_other1=sect_agriculture
      sect_other1=sf*(sect_agriculture+sect_agriculture2+sect_agriculture3+sect_solvents+sect_waste)
      sect_energy1=sf*sect_energy
      sect_industry1=sf*(sect_industry+sect_industry2)
      sect_transport1=sf*(sect_transport+sect_transport2)
      sect_resident1=sf*sect_resident
      sect_shipping1=sf*sect_shipping
      
      sect_other1= sect_other1* 1.e3 * avog/mw_s* 1.e-4
      sect_energy1 = sect_energy1 * 1.e3 * avog/mw_s* 1.e-4
      sect_industry1 = sect_industry1* 1.e3 * avog/mw_s* 1.e-4
      sect_transport1 = sect_transport1 * 1.e3 * avog/mw_s* 1.e-4
      sect_resident1 = sect_resident1 * 1.e3 * avog/mw_s* 1.e-4
      sect_shipping1 = sect_shipping1 * 1.e3 * avog/mw_s* 1.e-4

   else
      ; apply unit conversion for each species in the list separately assuming mw in file is correct
      ; since not correct, read correct mw:

      ; mw_s = mw_spec
      mw_s = mw_spec_set(ced)

      sect_other1_ced=sect_agriculture
      sect_other1_ced=sf*(sect_agriculture+sect_agriculture2+sect_agriculture3+sect_solvents+sect_waste)
      sect_energy1_ced=sect_energy
      sect_energy1_ced=sf*sect_energy
      sect_industry1_ced=sect_industry
      sect_industry1_ced=sf*(sect_industry+sect_industry2)
      sect_transport1_ced=sect_transport
      sect_transport1_ced=sf*(sect_transport+sect_transport2)
      sect_resident1_ced=sect_resident
      sect_resident1_ced=sf*sect_resident
      sect_shipping1_ced=sect_shipping
      sect_shipping1_ced=sf*sect_shipping
      
      sect_other1_ced= sect_other1_ced* 1.e3 * avog/mw_s* 1.e-4
      sect_energy1_ced = sect_energy1_ced * 1.e3 * avog/mw_s* 1.e-4
      sect_industry1_ced = sect_industry1_ced* 1.e3 * avog/mw_s* 1.e-4
      sect_transport1_ced = sect_transport1_ced * 1.e3 * avog/mw_s* 1.e-4
      sect_resident1_ced = sect_resident1_ced * 1.e3 * avog/mw_s* 1.e-4
      sect_shipping1_ced = sect_shipping1_ced * 1.e3 * avog/mw_s* 1.e-4


      sect_other1=sect_other1 + sect_other1_ced 
      sect_energy1=sect_energy1 +  sect_energy1_ced
      sect_industry1=sect_industry1 + sect_industry1_ced
      sect_transport1=sect_transport1 + sect_transport1_ced
      sect_resident1=sect_resident1 + sect_resident1_ced
      sect_shipping1=sect_shipping1 + sect_shipping1_ced
   end if  
     
 delete(files)
end do ;ced

      other_out1=area_conserve_remap_Wrap(lon_inp,lat_inp,sect_other1,lon_out,lat_out,False)
      energy_out1=area_conserve_remap_Wrap(lon_inp,lat_inp,sect_energy1,lon_out,lat_out,False)
      industry_out1=area_conserve_remap_Wrap(lon_inp,lat_inp,sect_industry1,lon_out,lat_out,False)
      transport_out1=area_conserve_remap_Wrap(lon_inp,lat_inp,sect_transport1,lon_out,lat_out,False)
      resident_out1=area_conserve_remap_Wrap(lon_inp,lat_inp,sect_resident1,lon_out,lat_out,False)
      shipping_out1=area_conserve_remap_Wrap(lon_inp,lat_inp,sect_shipping1,lon_out,lat_out,False)

print(min(energy_out1))
print(max(energy_out1))
print(min(industry_out1))
print(max(industry_out1))


;     sector = (/"other","energy","industry","transport","resident","shipping"/)

      other_out  = new((/ntim,nlat_out,nlon_out/),float)
      other_out!0 = "time"
      other_out!1 = "lat"
      other_out!2 = "lon"
      other_out&time = time
      other_out&lat= lat_out 
      other_out&lon= lon_out
      other_out = (/doubletofloat(other_out1)/)
      other_out@units = "molecules/cm2/s"
      other_out@molecular_weight=mw
      other_out@molecular_weight_units = "g/mole"
      other_out@short_name=other_out1@short_name 
      other_out@sector="agriculture,solvents,waste"


      energy_out  = new((/ntim,nlat_out,nlon_out/),float)
      energy_out!0 = "time"
      energy_out!1 = "lat"
      energy_out!2 = "lon"
      energy_out&time = time
      energy_out&lat= lat_out 
      energy_out&lon= lon_out
      energy_out@units = "molecules/cm2/s"
      energy_out@molecular_weight=mw
      energy_out@molecular_weight_units = "g/mole"
      energy_out@short_name=energy_out1@short_name 
      energy_out@sector="energy"
      energy_out = (/doubletofloat(energy_out1)/)

      industry_out  = new((/ntim,nlat_out,nlon_out/),float)
      industry_out!0 = "time"
      industry_out!1 = "lat"
      industry_out!2 = "lon"
      industry_out&time = time
      industry_out&lat= lat_out 
      industry_out&lon= lon_out
      industry_out@units = "molecules/cm2/s"
      industry_out@molecular_weight=mw
      industry_out@molecular_weight_units = "g/mole"
      industry_out@short_name=industry_out1@short_name 
      industry_out@sector="industry"
      industry_out = (/doubletofloat(industry_out1)/)

      resident_out  = new((/ntim,nlat_out,nlon_out/),float)
      resident_out!0 = "time"
      resident_out!1 = "lat"
      resident_out!2 = "lon"
      resident_out&time = time
      resident_out&lat= lat_out 
      resident_out&lon= lon_out
      resident_out@units = "molecules/cm2/s"
      resident_out@molecular_weight=mw
      resident_out@molecular_weight_units = "g/mole"
      resident_out@short_name=resident_out1@short_name 
      resident_out@sector="resident"
      resident_out = (/doubletofloat(resident_out1)/)


      transport_out  = new((/ntim,nlat_out,nlon_out/),float)
      transport_out!0 = "time"
      transport_out!1 = "lat"
      transport_out!2 = "lon"
      transport_out&time = time
      transport_out&lat= lat_out 
      transport_out&lon= lon_out
      transport_out@units = "molecules/cm2/s"
      transport_out@molecular_weight=mw
      transport_out@molecular_weight_units = "g/mole"
      transport_out@short_name=transport_out1@short_name 
      transport_out@sector="transport"
      transport_out = (/doubletofloat(transport_out1)/)

      shipping_out  = new((/ntim,nlat_out,nlon_out/),float)
      shipping_out!0 = "time"
      shipping_out!1 = "lat"
      shipping_out!2 = "lon"
      shipping_out&time = time
      shipping_out&lat= lat_out 
      shipping_out&lon= lon_out
      shipping_out@units = "molecules/cm2/s"
      shipping_out@molecular_weight=mw
      shipping_out@molecular_weight_units = "g/mole"
      shipping_out@short_name=shipping_out1@short_name 
      shipping_out@sector="shipping"
      shipping_out = (/doubletofloat(shipping_out1)/)

      

      ; output to netcdf file
;
if nspc .le. 5  then
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+spc_name2 + "_"+version+ "_c" + tostring(creation_date) + ".nc"
else
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+spc_name + "_"+version+ "_c" + tostring(creation_date) + ".nc"
end if
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

printVarSummary(other_out)
print(min(other_out))
print(max(other_out))
print(min(energy_out))
print(max(energy_out))
      ncdf->emiss_other= other_out
      ncdf->emiss_ene= energy_out
      ncdf->emiss_ind= industry_out
      ncdf->emiss_res= resident_out
      ncdf->emiss_ship= shipping_out
      ncdf->emiss_tran= transport_out
      ncdf->date         = date
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)
;;;;;;;;;;;

      
    if spc_name .eq. "BC" .or. spc_name .eq. "OC" then
        PI  = get_pi("float" )
   ; calculate number emissions
        diam = 0.134e-6
        if spc_name .eq. "BC" then
           rho = 1700.
        else
          rho = 1000.
        end if
        mw = 12.
        mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
        other_out = other_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        energy_out = energy_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        industry_out = industry_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        resident_out = resident_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        shipping_out = shipping_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        transport_out = transport_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
   
      other_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"
      energy_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"
      industry_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"
      shipping_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"
      transport_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"
      resident_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"

        if spc_name .eq. "BC" then
         other_out@long_name =  "number of BC Emissions"
         energy_out@long_name =  "number of BC Emissions"
         industry_out@long_name =  "number of BC Emissions"
         shipping_out@long_name =  "number of BC Emissions"
         transport_out@long_name =  "number of BC Emissions"
         resident_out@long_name =  "number of BC Emissions"
        end if
        if spc_name .eq. "OC" then
         other_out@long_name =  "number of (1.4*OC) POM Emissions"
         energy_out@long_name =  "number of (1.4*OC) POM Emissions"
         industry_out@long_name =  "number of (1.4*OC) POM Emissions"
         shipping_out@long_name =  "number of (1.4*OC) POM Emissions"
         transport_out@long_name =  "number of (1.4*OC) POM Emissions"
         resident_out@long_name =  "number of (1.4*OC) POM Emissions"
        end if

; output to netcdf file
;
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_num_"+spc_name2 + "_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of number "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->num_emiss_other= other_out
      ncdf->num_emiss_ene= energy_out
      ncdf->num_emiss_ind= industry_out
      ncdf->num_emiss_res= resident_out
      ncdf->num_emiss_ship= shipping_out
      ncdf->num_emiss_tran= transport_out
      ncdf->date         = date
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)
;;;;;;;;;;;

   end if ; bc and pom 


     if spc_name .eq. "SO2" then
      mw_so2 = 64.
      mw_so4_bam = 96.
      mw_so4_mam = 115.
      other_out_so2 = other_out
      other_out_so4 = other_out
      other_out_so2 = 0.975*other_out_so2
      other_out_so4 = 0.025*other_out_so4
      energy_out_so2 = energy_out
      energy_out_so4 = energy_out
      energy_out_so2 = 0.975*energy_out_so2
      energy_out_so4 = 0.025*energy_out_so4
      industry_out_so2 = industry_out
      industry_out_so4 = industry_out
      industry_out_so2 = 0.975*industry_out_so2
      industry_out_so4 = 0.025*industry_out_so4
      resident_out_so2 = resident_out 
      resident_out_so4 = resident_out 
      resident_out_so2 = 0.975*resident_out_so2
      resident_out_so4 = 0.025*resident_out_so4
      shipping_out_so2 = shipping_out
      shipping_out_so4 = shipping_out
      shipping_out_so2 = 0.975*shipping_out_so2
      shipping_out_so4 = 0.025*shipping_out_so4
      transport_out_so2 = transport_out
      transport_out_so4 = transport_out
      transport_out_so2 = 0.975*transport_out_so2
      transport_out_so4 = 0.025*transport_out_so4

      other_out_so2@long_name="0.975*SO2 Other Emissions" ; surface
      other_out_so4@long_name="0.025*SO2 Other Emissions" ; surface
      industry_out_so2@long_name="0.975*SO2 Industry Emissions" ; surface
      industry_out_so4@long_name="0.025*SO2 Industry Emissions" ; surface
      energy_out_so2@long_name="0.975*SO2 Energy Emissions" ; surface
      energy_out_so4@long_name="0.025*SO2 Energy Emissions" ; surface
      resident_out_so2@long_name="0.975*SO2 Resident Emissions" ; surface
      resident_out_so4@long_name="0.025*SO2 Resident Emissions" ; surface
      shipping_out_so2@long_name="0.975*SO2 Shipping Emissions" ; surface
      shipping_out_so4@long_name="0.025*SO2 Shipping Emissions" ; surface
      transport_out_so2@long_name="0.975*SO2 Transport Emissions" ; surface
      transport_out_so4@long_name="0.025*SO2 Transport Emissions" ; surface
      
       nlon = dimsizes(lon_out)
       nlat = dimsizes(lat_out)
       ; write so4_a1: sulfate accumulation mode for energy+industrial
       ; vertically distributed: 100-300m
       nalt = 8
       dz = 0.05 ;km
       dalt = 2.e4 ;cm = 200m
       altitude = fspan(0,nalt-1,nalt) * dz + 0.5*dz
       altitude_int = fspan(0,nalt,nalt+1) * dz
       altitude!0="altitude"
       altitude_int!0="altitude_int"
       altitude@units="km"
       altitude@long_name="Altitude"
       altitude_int@units="km"
       altitude_int@long_name="Altitude interfaces"
       energy_out_so4_vert = new((/ntim,nalt,nlat,nlon/),"float")
       industry_out_so4_vert = new((/ntim,nalt,nlat,nlon/),"float")
       energy_out_so4_vert(:,:,:,:) = 0.
       industry_out_so4_vert(:,:,:,:) = 0.
       energy_out_so4_vol = energy_out_so4
       energy_out_so4_vol = energy_out_so4/dalt
       industry_out_so4_vol = industry_out_so4
       industry_out_so4_vol = industry_out_so4/dalt
       do ialt=2,5
         energy_out_so4_vert(:,ialt,:,:) = energy_out_so4_vol
         industry_out_so4_vert(:,ialt,:,:) = industry_out_so4_vol
       end do
       energy_out_so4_vert!0="time"
       energy_out_so4_vert!1="altitude"
       energy_out_so4_vert!2="lat"
       energy_out_so4_vert!3="lon"
       energy_out_so4_vert@units="kg m-3 s-1"
       energy_out_so4_vert@long_name="0.025*SO2 Anthropogenic Emissions from Energy"
       energy_out_so4_vert&altitude=altitude

       industry_out_so4_vert!0="time"
       industry_out_so4_vert!1="altitude"
       industry_out_so4_vert!2="lat"
       industry_out_so4_vert!3="lon"
       industry_out_so4_vert@units="kg m-3 s-1"
       industry_out_so4_vert@long_name="0.025*SO2 Anthropogenic Emissions from Industry"
       industry_out_so4_vert&altitude=altitude

printVarSummary(energy_out_so4_vert)
printVarSummary(industry_out_so4_vert)

; spec_name_sec = (/"SO2_anthro-ag-ship-res_surface","SO2_anthro-ene_surface","so4_a1_anthro-ene_vertical","so4_a1_anthro-ag-ship_surface","so4_a2_anthro-res_surface","num_so4_a1_anthro-ene_vertical","num_so4_a1_anthro-ag-ship_surface","num_so4_a2_anthro-res_surface"/)
      printVarSummary(time)
      

  
   ; calculate number emission for so4 
   ;  sector = (/"other","energy","industry","transport","resident","shipping"/)
     ;varname = 'emiss_ag_sol_was'
        diam = 0.134e-6
        rho = 1770.
        mw2 = 115.
        mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
        num_other_out_so4 = other_out_so4 
        num_other_out_so4= num_other_out_so4*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_other_out_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_other_out_so4@long_name = "particle number emissions of SO2*0.025 "

    ; if sec.eq.1 .or. sec.eq.2 .or.  sec.eq.5 then
    ; varname = 'emiss_shipping'
    ;varname = 'emiss_ene_ind'
        diam = 0.261e-6
        rho = 1770.
        mw2 = 115.
        num_energy_out_so4_vert = energy_out_so4_vert
        num_energy_out_so4_vert= num_energy_out_so4_vert*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_energy_out_so4_vert@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_energy_out_so4_vert@long_name = "particle number emissions of SO2*0.025 "
        num_industry_out_so4_vert = industry_out_so4_vert
        num_industry_out_so4_vert= num_industry_out_so4_vert*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_industry_out_so4_vert@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_industry_out_so4_vert@long_name = "particle number emissions of SO2*0.025 "
        num_shipping_out_so4 = shipping_out_so4 
        num_shipping_out_so4= num_shipping_out_so4*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_shipping_out_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_shipping_out_so4@long_name = "particle number emissions of SO2*0.025 "
     ;end if
     ; if sec.eq.3 .or. sec.eq.4 then
     ;varname = 'emiss_res_tran'
        diam = 0.0504e-6
        rho = 1770.
        mw2 = 115.
        num_resident_out_so4= resident_out_so4
        num_resident_out_so4= num_resident_out_so4*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_resident_out_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_resident_out_so4@long_name = "particle number emissions of SO2*0.025 "
        num_transport_out_so4= transport_out_so4
        num_transport_out_so4= num_transport_out_so4*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_transport_out_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_transport_out_so4@long_name = "particle number emissions of SO2*0.025 "
     ; end if


      ; output to netcdf file
; spec_name_sec = (/"SO2_anthro-ag-ship-res_surface","SO2_anthro-ene_surface","so4_a1_anthro-ene_vertical","so4_a1_anthro-ag-ship_surface","so4_a2_anthro-res_surface","num_so4_a1_anthro-ene_vertical","num_so4_a1_anthro-ag-ship_surface","num_so4_a2_anthro-res_surface"/)
;
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_SO2_anthro_ag-ship-res-surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so2
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_other= other_out_so2
      ncdf->emiss_res= resident_out_so2
      ncdf->emiss_tran= transport_out_so2
      ncdf->emiss_ship= shipping_out_so2
      ncdf->date         = date
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_SO2_anthro_ene-surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so2
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_ene= energy_out_so2
      ncdf->emiss_ind= industry_out_so2
      ncdf->date         = date
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;
           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_so4_a1_anthro-ag-ship_surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_other= other_out_so4
      ncdf->emiss_ship= shipping_out_so4
      ncdf->date         = date
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;
           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_num_so4_a1_anthro-ag-ship_surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_other= num_other_out_so4
      ncdf->emiss_ship= num_shipping_out_so4
      ncdf->date         = date
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;
           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_so4_a2_anthro-res_surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_tran= transport_out_so4
      ncdf->emiss_res= resident_out_so4
      ncdf->date         = date
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_num_so4_a2_anthro-res_surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_tran= num_transport_out_so4
      ncdf->emiss_res= num_resident_out_so4
      ncdf->date         = date
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;



           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_so4_a1_anthro-ene_vertical_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_ene= energy_out_so4_vert
      ncdf->emiss_ind= industry_out_so4_vert
      ncdf->altitude = altitude
      ncdf->altitude_int = altitude_int
      ncdf->date         = date
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;

           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_num_so4_a1_anthro-ene_vertical_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)

      ncdf->emiss_ene= energy_out_so4_vert
      ncdf->emiss_ind= industry_out_so4_vert
      ncdf->altitude = altitude
      ncdf->altitude_int = altitude_int
      ncdf->date         = date
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;

           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_num_so4_a1_anthro-ene_vertical_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_ene= num_energy_out_so4_vert
      ncdf->emiss_ind= num_industry_out_so4_vert
      ncdf->date         = date
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      ncdf->altitude = altitude
      ncdf->altitude_int = altitude_int
      delete(ncdf)

;;;;;;;;;;;


   end if  ; SO2
;

;



;
; clean up
;
print("clean up")
      ;delete(wk)
      delete(in)
      delete(file_inp)
      delete(time)
      ;delete(time_bnds)
      delete(utc_date)
      delete(date)
;    delete(var_name)
    delete(mw)
    delete(spc_name)
; 
   ; end loop over sectors 
; end loop over species
;
else
do sec=0,0
 if spc_name.eq. "IVOC"
   specs_in = (/"C3H6", "C3H8", "C2H6", "C2H4", "BIGENE", "BIGALK", "CH3COCH3", "MEK", "CH3CHO", "CH2O", "BENZENE", "TOLUENE", "XYLENES"/)
   hist_spec = "IVOC=0.2*("+str_join(specs_in,"+")+")"
   fact = 0.2
 end if
  if spc_name.eq. "SVOC"
   specs_in = "pom_a4"
   hist_spec = "SVOC=0.6*pom_a4"
   fact = 0.6
  end if
 if spc_name.eq. "SOAGx1.5"
; if sec.eq.0 then
;  specs_in  = (/"BIGALK","BIGENE","TOLUENE","BENZENE","XYLENES"/)
;  yield  = (/0.05, 0.05, 0.15, 0.15, 0.15/)
;  hist_spec = "SOAG=1.5* yields(0.05, 0.05, 0.15, 0.15, 0.15) from ("+str_join(specs_in,"+")+")"
; end if
; if sec.eq.1 then
   specs_in  = (/"BIGALK","BIGENE","TOLUENE","BENZENE","XYLENES","ISOP","MTERP"/)
   yield  = (/0.05, 0.05, 0.15, 0.15, 0.15, 0.04, 0.25/)
   hist_spec = "SOAG=1.5* yields(0.05, 0.05, 0.15, 0.15, 0.15, 0.04, 0.25) from ("+str_join(specs_in,"+")+")"
; end if
   fact = 1.5
 end if
 ; read required files
   cnt_s = 0
   nsp = dimsizes(specs_in)
   do sp = 0,nsp-1
    files = systemfunc("ls " + dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+specs_in(sp)+"_"+version+ "_c*.nc")
     print(files)
    files1=files(0)
    in = addfile(files1,"r")

    if cnt_s.eq.0 then
     time1 = in->time
      date = in->date
      lat = in->lat
      lon = in->lon
      emiss_o = in->emiss_other
      emiss_e = in->emiss_ene
      emiss_i = in->emiss_ind
      emiss_t = in->emiss_tran
      emiss_r = in->emiss_res
      emiss_s = in->emiss_ship

      mw_hc =  emiss_o@molecular_weight
     if spc_name.eq."SOAGx1.5" then
      emiss_o = yield(sp)*emiss_o*mw_hc/mw
      emiss_e = yield(sp)*emiss_e*mw_hc/mw
      emiss_i = yield(sp)*emiss_i*mw_hc/mw
      emiss_t = yield(sp)*emiss_t*mw_hc/mw
      emiss_r = yield(sp)*emiss_r*mw_hc/mw
      emiss_s = yield(sp)*emiss_s*mw_hc/mw
     else
      emiss_o = emiss_o *mw_hc/mw
      emiss_e = emiss_e *mw_hc/mw
      emiss_i = emiss_i *mw_hc/mw
      emiss_t = emiss_t *mw_hc/mw
      emiss_r = emiss_r *mw_hc/mw
      emiss_s = emiss_s *mw_hc/mw
     end if
    else
      emiss_o_hc1 = in->emiss_other
      emiss_e_hc1 = in->emiss_ene
      emiss_i_hc1 = in->emiss_ind
      emiss_t_hc1 = in->emiss_tran
      emiss_r_hc1 = in->emiss_res
      emiss_s_hc1 = in->emiss_ship
      mw_hc =  emiss_o_hc1@molecular_weight
     if spc_name.eq."SOAGx1.5" then
       emiss_o = emiss_o + yield(sp)*emiss_o_hc1*mw_hc/mw
       emiss_e = emiss_e + yield(sp)*emiss_e_hc1*mw_hc/mw
       emiss_i = emiss_i + yield(sp)*emiss_i_hc1*mw_hc/mw
       emiss_t = emiss_t + yield(sp)*emiss_t_hc1*mw_hc/mw
       emiss_r = emiss_r + yield(sp)*emiss_r_hc1*mw_hc/mw
       emiss_s = emiss_s + yield(sp)*emiss_s_hc1*mw_hc/mw
     else
       emiss_o = emiss_o + emiss_o_hc1*mw_hc/mw
       emiss_e = emiss_e + emiss_e_hc1*mw_hc/mw
       emiss_i = emiss_i + emiss_i_hc1*mw_hc/mw
       emiss_t = emiss_t + emiss_t_hc1*mw_hc/mw
       emiss_r = emiss_r + emiss_r_hc1*mw_hc/mw
       emiss_s = emiss_s + emiss_s_hc1*mw_hc/mw
     end if
      delete(emiss_o_hc1)
      delete(emiss_e_hc1)
      delete(emiss_i_hc1)
      delete(emiss_t_hc1)
      delete(emiss_r_hc1)
      delete(emiss_s_hc1)
    end if
    cnt_s = cnt_s+1
   end do
    ;IVOC = 0.2*(HCs)
    ;SVOC = 0.6*pom_a4
     emiss_o = fact*emiss_o
     emiss_e = fact*emiss_e
     emiss_i = fact*emiss_i
     emiss_t = fact*emiss_t
     emiss_r = fact*emiss_r
     emiss_s = fact*emiss_s

      emiss_other = emiss_o
      emiss_other@history= hist_spec
      emiss_other@long_name= spc_name+" CAMS emissions"
      emiss_ene = emiss_e 
      emiss_ene@history= hist_spec
      emiss_ene@long_name= spc_name+" CAMS emissions"
      emiss_ind = emiss_i
      emiss_ind@history= hist_spec
      emiss_ind@long_name= spc_name+" CAMS emissions"
      emiss_tran = emiss_t
      emiss_tran@history= hist_spec
      emiss_tran@long_name= spc_name+" CAMS emissions"
      emiss_res = emiss_r
      emiss_res@history= hist_spec
      emiss_res@long_name= spc_name+" CAMS emissions"
      emiss_ship  = emiss_s
      emiss_ship@history= hist_spec
      emiss_ship@long_name= spc_name+" CAMS emissions"

  file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+spc_name + version+ "_c" + tostring(creation_date) + ".nc"
  system("/bin/rm -f " + file_out) ; remove any pre-existing file
  ncdf = addfile(file_out,"c")
  system("/bin/rm -f " + file_out) ; remove any pre-existing file
  ncdf = addfile(file_out,"c")

      file_att               = True
      file_att@molecular_weight = mw
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@creation_date = creation_date
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_se.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
      ncdf->emiss_other= emiss_other
      ncdf->emiss_ene= emiss_ene
      ncdf->emiss_ind= emiss_ind
      ncdf->emiss_res= emiss_res
      ncdf->emiss_ship= emiss_ship
      ncdf->emiss_tran= emiss_tran

  ncdf->date         = date
  ncdf->time         = time1
  ncdf->lat          = lat
  ncdf->lon          = lon

  delete(ncdf)
  delete(time1)
  delete(date)
  delete(specs_in)

 if spc_name.eq."SOAGx1.5" then
  delete(yield)
 end if
 end do; SOAG, IVOC, SVOC

end if
  delete(spc_name1)
 end do  ;  species
;
end
