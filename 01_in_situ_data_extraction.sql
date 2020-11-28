/*
VARIABLES TO EXTRACT:

cruise id (mission name)
mission descriptor
event id
sample_id
station
latitude
longitude
year
month
day
time
depth
parameter name
value
method
quality code


PARAMETERS TO EXTRACT:
(source: 01_BC_Data_Retrievals_08Mar2017.txt)

90000160 | 90 | HPLC_ACAROT | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | alpha-carotene
90000161 | 90 | HPLC_ALLOX | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | alloxanthin
90000162 | 90 | HPLC_ASTAX | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | astaxanthin
90000163 | 90 | HPLC_BCAROT | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | beta-carotene
90000165 | 90 | HPLC_BUTLIKE | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | hplc pigment similar to but19
90000164 | 90 | HPLC_BUT19 | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | 19-butanoyl-oxy-fucoxanthin
90000166 | 90 | HPLC_CHLA | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | hplc chlorophyll a (includes chlorophyllide a+)
90000167 | 90 | HPLC_CHLB | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | chlorophyll-b (includes dv-chlb if any present)
90000168 | 90 | HPLC_CHLC12 | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | chlorophyll - c1+c2
90000182 | 90 | HPLC_CHLC3 | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | chlorophyll - c3
90000169 | 90 | HPLC_CHLIDEA | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | chlorophyllide-a
90000170 | 90 | HPLC_DIADINOX | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | diadinoxanthin
90000171 | 90 | HPLC_DIATOX | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | diatoxanthin
90000289 | 20 | HPLC_%DVa | % | 3 | 3 | 0.00000 | 100.00000 | MARY | Percent DVa
90000290 | 20 | HPLC_%DVb | % | 3 | 3 | 0.00000 | 100.00000 | MARY | Percent DVb
90000172 | 90 | HPLC_FUCOX | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | fucoxanthin
90000174 | 90 | HPLC_HEXLIKE | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | hplc pigment similar to hex19 (elutes ~0.1 minute after 19hex)
90000181 | 90 | HPLC_HEXLIKE2 | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | hplc pigment similar to hex19 (elutes ~0.2 minute after 19hex)
90000173 | 90 | HPLC_HEX19 | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | 19-hexanoyl-oxy-fucoxanthin
90000175 | 90 | HPLC_PERID | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | peridinin
90000176 | 90 | HPLC_PHAEO | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | total hplc phaeopigment concentration (includes all phaeophorbides and phaeophytins)
90000180 | 90 | HPLC_PRASINOX | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | prasinoxanthin
90000177 | 90 | HPLC_PYROPHAE | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | pytophaeophorbide a
90000178 | 90 | HPLC_VIOLAX | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | violaxanthin
90000179 | 90 | HPLC_ZEA | mg/m3 | 2 | 3 | 0.00000 | 64.00000 | KennedyM | zeaxanthin

90000003 | 90 | Temperature | degrees C | 3 | 3 | -3.00000 | 35.00000 | P. Strain | Temperature
90000029 | 90 | Chlorophyll A | mg/m3 | 3 | 6 | 0.00000 | 200.00000 | P. Strain | Chlorophyll a
90000002 | 90 | Salinity | none | 3 | 3 | 0.00000 | 400.00000 | P. Strain | Salinity
90000084 | 90 | PAR | uE m-2 s-1 | 4 | 3 | 0.00000 | 5000.00000 | P. Strain | Photosynthetically active radiation (400-700 nm)
90000014 | 90 | Nitrate | mmol/m3 | 4 | 2 | -1.00000 | 75.00000 | P. Strain | Nitrate + nitrite
90000013 | 90 | Phosphate | mmol/m3 | 4 | 3 | 0.00000 | 50.00000 | P. Strain | Dissolved inorganic phosphate
90000079 | 90 | Phosphate: part. | mmol/m3 | 3 | 3 | 0.00000 | 30.00000 | P. Strain | Particulate phosphate
90000012 | 90 | Silicate | mmol/m3 | 4 | 2 | 0.00000 | 260.00000 | P. Strain | Dissolved reactive silicate
*/


SELECT DISTINCT
   mi.name cruise_id,
   mi.descriptor mission_descriptor,
   ev.collector_event_id event_id,
   dh.collector_sample_id sample_id,
   ev.collector_station_name station,
   dh.start_lat latitude,
   dh.start_lon longitude,
   TO_NUMBER(TO_CHAR(dh.start_date, 'YYYY')) year,
   TO_NUMBER(TO_CHAR(dh.start_date, 'MM')) month,
   TO_NUMBER(TO_CHAR(dh.start_date, 'DD')) day,
   dh.start_time time,
   dh.start_depth,
   dr.parameter_name parameter,
   dd.data_value value,
   dt.method,
   dd.data_qc_code QC
FROM
   biochem.bcdiscretedtails dd,
   biochem.bcdiscretehedrs dh,
   biochem.bcdataretrievals dr,
   biochem.bcdatatypes dt,
   biochem.bcevents ev,
   biochem.bcmissions mi
WHERE
   /* link tables - only include results from parent tables that are connected to records in child tables */
   mi.mission_seq = ev.mission_seq
   AND ev.event_seq = dh.event_seq
   AND dh.discrete_seq = dd.discrete_seq
   AND dr.data_retrieval_seq = dt.data_retrieval_seq
   AND dt.data_type_seq = dd.data_type_seq
   /* non null data only */
   AND dd.data_value IS NOT NULL
   
   /* optional year, latitude, and longitude bounds, with formatting examples - note: BETWEEN is inclusive */
   
   AND TO_NUMBER(TO_CHAR(dh.start_date, 'YYYY')) BETWEEN 1997 AND 2020
   /* AND TO_NUMBER(TO_CHAR(dh.start_date, 'YYYY')) >= 1995 */
   /* AND dh.start_lat BETWEEN 39 AND 82 */
   /* AND dh.start_lat >= 39 */
   /* AND dh.start_lon BETWEEN -95 AND -41 */
   /* AND dh.start_lon >= -95 */
   
   /* data type filter for HPLC, from BC_Data_Retrievals_08Mar2017.txt */
   AND dt.data_retrieval_seq IN ('90000160','90000161','90000162','90000163','90000165','90000164','90000166','90000167','90000168','90000182','90000169','90000170','90000171','90000289','90000290','90000172','90000174','90000181','90000173','90000175','90000176','90000180','90000177','90000178','90000179','90000003','90000029','90000002','90000084','90000014','90000013','90000079','90000012')
ORDER BY
   year,
   month,
   day,
   dh.start_depth
