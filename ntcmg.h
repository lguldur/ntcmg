/*
ntcmg V1.0

This code implements the NTCM-G algorithm as described
in the "EUROPEAN GNSS (GALILEO) OPEN SERVICE/NTCM-G IONOSPHERIC MODEL DESCRIPTION" document:
https://www.gsc-europa.eu/sites/default/files/NTCM-G_Ionospheric_Model_Description_-_v1.0.pdf

Copyright (C) 2023 by David Duchet
This code is licensed under MIT license (see LICENSE.txt for details)

(dd@bruine.org, FRSEC18E62CA1787)
*/


namespace ntcmg
{
double stec( // TECU
	double ai0,
	double ai1,
	double ai2,
	double rxlat, // rad
	double rxlon, // rad
	double rxalt, // m
	double satlat, // rad
	double satlon, // rad
	double satalt,// m
	double utc_time, // hour
	double doy // day
	);


double iono_delay(
	double STEC, // TECU
	double freq // Hz
	);


double iono_delay( // TECU
	double ai0,
	double ai1,
	double ai2,
	double freq, // Hz
	double rxlat, // rad
	double rxlon, // rad
	double rxalt, // m
	double satlat, // rad
	double satlon, // rad
	double satalt,// m
	double utc_time, // hour
	double doy // day
	);


double test(); // 0 is good. ~0.0003 is  good enough ;-)

};// namespace
