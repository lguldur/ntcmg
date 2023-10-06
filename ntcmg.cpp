/*
ntcmg V1.0

This code implements the NTCM-G algorithm as described
in the "EUROPEAN GNSS (GALILEO) OPEN SERVICE/NTCM-G IONOSPHERIC MODEL DESCRIPTION" document:
https://www.gsc-europa.eu/sites/default/files/NTCM-G_Ionospheric_Model_Description_-_v1.0.pdf

Copyright (C) 2023 by David Duchet
This code is licensed under MIT license (see LICENSE file for details)

(dd@bruine.org, FRSEC18E62CA1787)
*/


#include <cmath>

namespace ntcmg
{

// constants
double a=6378137.; // m
double b=6356752.3142; // m
double pi=3.1415926535898; // rad
double Re=6371.; // km
double hI=450.; // km
double PHIgnp=79.74; // deg
double LAMBDAgnp=-71.78; // deg


void lla2xyz(
	double lat, // rad
	double lon, // rad
	double alt, // rad
	double & x, // m
	double & y, // m
	double & z // m
	)
{
	double e=sqrt(1.-b*b/(a*a));
	double e2=e*e;
	
	double cosphi=cos(lat);
	double sinphi=sin(lat);
	double coslambda=cos(lon);
	double sinlambda=sin(lon);
	double Nphi=a/sqrt(1-e2*sinphi*sinphi); // Eq. 19
	
	x=(Nphi+alt)*cosphi*coslambda; // Eq. 16
	y=(Nphi+alt)*cosphi*sinlambda; // Eq. 17
	z=(Nphi*(1-e2)+alt)*sinphi; // Eq. 18
}


void elevazim(
	double rxx, // m
	double rxy, // m
	double rxz, // m
	double satx, // m
	double saty, // m
	double satz, // m
	double rxlat, // rad
	double rxlon, // rad
	double & elev, // rad
	double & azim // rad	
	)
{
	// Eq. 20
	double dx=satx-rxx;
	double dy=saty-rxy;
	double dz=satz-rxz;
	
	// Eq. 21
	double cosphi=cos(rxlat);
	double sinphi=sin(rxlat);
	double coslambda=cos(rxlon);
	double sinlambda=sin(rxlon);
	
	double m11=-sinphi*coslambda; 	double m12=-sinphi*sinlambda;	double m13=cosphi;
	double m21=-sinlambda;			double m22=coslambda;			double m23=0;
	double m31=cosphi*coslambda;	double m32=cosphi*sinlambda;	double m33=sinphi;
	
	double Dx=m11*dx + m12*dy + m13*dz;
	double Dy=m21*dx + m22*dy + m23*dz;
	double Dz=m31*dx + m32*dy + m33*dz;

	// Eq. 22
	azim=atan2(Dy, Dx);
	if (azim < 0)
		azim+=2.*pi;
	
	elev = (0.5 - atan2(sqrt(Dx * Dx + Dy * Dy), Dz) / pi) * pi; // Eq. 23
}


void ipp(
	double elev, // rad
	double azim, // rad
	double rxlat, // rad
	double rxlon, // rad
	double & ipplat, // rad
	double & ipplon // rad
	)
{
	double psipp=pi/2.-elev-asin(Re*cos(elev)/(Re+hI)); // Eq. 24
	ipplat=asin(sin(rxlat)*cos(psipp)+cos(rxlat)*sin(psipp)*cos(azim)); // Eq. 25
	ipplon=rxlon+asin(sin(psipp)*sin(azim)/cos(ipplat)); // Eq. 26
}


double localtime( // hours
	double ipplon, // rad
	double utc_time // hour
	)
{
	return utc_time+ipplon*180./pi/15.; // Eq. 27
}


double sun_declination( // rad
	double doy // day
	)
{
	return 23.44 * sin(0.9856* (doy-80.7)*pi/180.)*pi/180.; // Eq. 28
}


void solar_zenith_angle_dependence(
	double ipplat, // rad
	double sun_declination, // rad
	double & coskhi2, // -
	double & coskhi3 // -
	)
{
	double Pf1=0.4;
	
	coskhi3=cos(ipplat-sun_declination)+Pf1; // Eq. 29
	coskhi2=cos(ipplat-sun_declination)-2.*ipplat*sin(sun_declination)/pi; // Eq. 30
}


double geomagnetic_ipplat( // rad
	double ipplat, // rad
	double ipplon // rad
)
{
	// Eq. 31
	double ipplatsin = sin(ipplat) * sin(PHIgnp*pi/180.) + cos(ipplat) * cos(PHIgnp * pi / 180.) * cos(ipplon - LAMBDAgnp * pi / 180.);
	return asin(ipplatsin);
}


double vtec2stec_mapping_function(
	double elev // rad
	)
{
	double sinz=Re*sin(0.9782*(pi/2.-elev))/(Re+hI); // Eq. 33
	double MF=1./sqrt(1.-sinz*sinz); // Eq. 32
	return MF;
}


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
	)
{
	// iono k parameters
	double k1=0.92519;
	double k2=0.16951;
	double k3=0.00443;
	double k4=0.06626;
	double k5=0.00899;
	double k6=0.21289;
	double k7=-0.15414;
	double k8=-0.38439;
	double k9=1.14023;
	double k10=1.20556;
	double k11=1.41808; // TECU
	double k12=0.13985; // TECU/sfu
	
	double LTd=14.; // hour
	double doya=18.; // day
	double doysa=6.; // day
	double latc1deg=16.; // deg
	double latc2deg=-10.; // deg
	double sigmac1deg=12.; // deg
	double sigmac2deg=13.; // deg
	
	double latc1 = latc1deg *pi / 180.;
	double latc2 = latc2deg *pi / 180.;
	double sigmac1 = sigmac1deg *pi / 180.;
	double sigmac2 = sigmac2deg *pi / 180.;

	double twopi=2*pi;

	double rxx;
	double rxy;
	double rxz;
	double satx;
	double saty;
	double satz;

	lla2xyz(rxlat, rxlon, rxalt, rxx, rxy, rxz);
	lla2xyz(satlat, satlon, satalt, satx, saty, satz);

	double elev;
	double azim;

	elevazim(rxx, rxy, rxz, satx, saty, satz, rxlat, rxlon, elev, azim);

	double ipplat;
	double ipplon;

	ipp(elev, azim, rxlat, rxlon, ipplat, ipplon);
	
	// F1
	double LT = localtime(ipplon, utc_time);
	double Vd=twopi*(LT-LTd)/24.; // Eq. 5
	double Vsd=twopi*LT/12.; // Eq. 6
	double Vtd=twopi*LT/8.; // Eq. 7
	double delta = sun_declination(doy);
	double coskhi2;
	double coskhi3;
	solar_zenith_angle_dependence(ipplat, delta, coskhi2, coskhi3);
	double F1=coskhi3 + coskhi2 * ( k1*cos(Vd) + k2*cos(Vsd) + k3*sin(Vsd) + k4*cos(Vtd) + k5*sin(Vtd)); // Eq. 4
	
	// F2
	double Va=twopi*(doy-doya)/365.25; // Eq. 9
	double Vsa=2*twopi*(doy-doysa)/365.25; // Eq. 10
	double F2=1 + k6*cos(Va) + k7*cos(Vsa); // Eq. 8
	
	// F3
	double mag_ipplat = geomagnetic_ipplat(ipplat, ipplon);
	double F3=1 + k8*cos(mag_ipplat); // Eq. 11
	
	// F4
	double dlat1 = mag_ipplat - latc1;
	double dlat2 = mag_ipplat - latc2;
	double EC1=-dlat1*dlat1/(2*sigmac1*sigmac1); // Eq. 13
	double EC2=-dlat2*dlat2/(2*sigmac2*sigmac2); // Eq. 14
	double F4=1 + k9*exp(EC1) + k10*exp(EC2); // Eq. 12
	
	// F5
	double Azpar=sqrt(ai0*ai0 + 1633.33*ai1*ai1 + 4802000*ai2*ai2 + 3266.67*ai0*ai2); // Eq. 2
	double F5=k11 + k12*Azpar; // Eq. 15
	
	// VTEC NTCM-G
	double VTEC= F1*F2*F3*F4*F5; // Eq. 3
	
	// STEC
	return vtec2stec_mapping_function(elev)*VTEC; // Eq. 34
}


double iono_delay( // m
	double STEC,	// TECU
	double freq // Hz
	)
{
	return 40.3*STEC/(freq*freq); // Eq. 1
}


double iono_delay( // m
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
	)
{
	return iono_delay(stec(ai0, ai1, ai2, rxlat, rxlon, rxalt, satlat, satlon, satalt, utc_time, doy), freq);
}


double test() // 0 is good
{
	double high[]=
	{
		105,
		0,
		-62.34,
		82.49,
		78.11,
		8.23,
		54.29,
		20281546.18,
		33.7567,
		105,
		0,
		-62.34,
		82.49,
		78.11,
		-158.03,
		24.05,
		20275295.43,
		98.7476,
		105,
		0,
		-62.34,
		82.49,
		78.11,
		-30.86,
		41.04,
		19953770.93,
		42.2768,
		105,
		4,
		-62.34,
		82.49,
		78.11,
		-85.72,
		53.69,
		20544786.65,
		32.1058,
		105,
		4,
		-62.34,
		82.49,
		78.11,
		-130.77,
		54.40,
		20121312.46,
		33.6307,
		105,
		4,
		-62.34,
		82.49,
		78.11,
		140.68,
		35.85,
		19953735.00,
		54.6962,
		105,
		8,
		-62.34,
		82.49,
		78.11,
		-126.28,
		51.26,
		20513440.10,
		30.8902,
		105,
		8,
		-62.34,
		82.49,
		78.11,
		84.26,
		54.68,
		20305726.79,
		33.8241,
		105,
		8,
		-62.34,
		82.49,
		78.11,
		-96.21,
		37.33,
		19956072.48,
		36.9695,
		105,
		12,
		-62.34,
		82.49,
		78.11,
		81.09,
		35.20,
		20278071.03,
		65.0500,
		105,
		12,
		-62.34,
		82.49,
		78.11,
		175.57,
		51.89,
		19995445.72,
		33.9561,
		105,
		12,
		-62.34,
		82.49,
		78.11,
		4.25,
		53.43,
		20107681.66,
		40.4107,
		105,
		16,
		-62.34,
		82.49,
		78.11,
		14.89,
		32.88,
		20636367.33,
		73.7953,
		105,
		16,
		-62.34,
		82.49,
		78.11,
		-70.26,
		50.63,
		20043030.82,
		47.1845,
		105,
		16,
		-62.34,
		82.49,
		78.11,
		-130.60,
		49.21,
		20288021.34,
		46.4704,
		105,
		20,
		-62.34,
		82.49,
		78.11,
		-52.46,
		24.28,
		19831557.96,
		87.9216,
		105,
		20,
		-62.34,
		82.49,
		78.11,
		-165.78,
		35.06,
		20196268.24,
		69.7195,
		105,
		20,
		-62.34,
		82.49,
		78.11,
		168.73,
		52.58,
		20288372.95,
		44.4995,
		105,
		0,
		-52.81,
		5.25,
		-25.76,
		-89.48,
		-29.05,
		20081457.33,
		125.2377,
		105,
		0,
		-52.81,
		5.25,
		-25.76,
		-46.73,
		-24.08,
		19975517.42,
		86.3239,
		105,
		0,
		-52.81,
		5.25,
		-25.76,
		-99.26,
		34.47,
		20275286.46,
		127.0015,
		105,
		4,
		-52.81,
		5.25,
		-25.76,
		-46.61,
		54.84,
		20258938.89,
		75.5439,
		105,
		4,
		-52.81,
		5.25,
		-25.76,
		-85.72,
		53.68,
		20544786.61,
		86.2081,
		105,
		4,
		-52.81,
		5.25,
		-25.76,
		-18.13,
		14.17,
		20267783.18,
		59.7016,
		105,
		8,
		-52.81,
		5.25,
		-25.76,
		7.14,
		-19.55,
		20226657.45,
		83.0316,
		105,
		8,
		-52.81,
		5.25,
		-25.76,
		-48.38,
		-31.04,
		20069586.93,
		39.9556,
		105,
		8,
		-52.81,
		5.25,
		-25.76,
		-58.59,
		21.93,
		20008556.82,
		31.4675,
		105,
		12,
		-52.81,
		5.25,
		-25.76,
		-102.83,
		-40.74,
		20153844.84,
		194.4159,
		105,
		12,
		-52.81,
		5.25,
		-25.76,
		-0.60,
		10.75,
		20272829.17,
		181.2492,
		105,
		12,
		-52.81,
		5.25,
		-25.76,
		-120.35,
		11.00,
		20283503.35,
		184.0749,
		105,
		16,
		-52.81,
		5.25,
		-25.76,
		-70.26,
		50.63,
		20043030.72,
		216.2278,
		105,
		16,
		-52.81,
		5.25,
		-25.76,
		-72.73,
		-9.78,
		19936049.27,
		162.4852,
		105,
		16,
		-52.81,
		5.25,
		-25.76,
		-66.77,
		2.37,
		19986966.89,
		147.7473,
		105,
		20,
		-52.81,
		5.25,
		-25.76,
		-1.57,
		-7.90,
		20373709.74,
		220.0012,
		105,
		20,
		-52.81,
		5.25,
		-25.76,
		0.44,
		50.83,
		19975412.45,
		240.8009,
		105,
		20,
		-52.81,
		5.25,
		-25.76,
		10.94,
		44.72,
		20450566.19,
		252.0204,
	};
	
	double med[]=
	{
		105,
		0,
		40.19,
		-3.00,
		-23.32,
		76.65,
		-41.43,
		20157673.93,
		28.3208,
		105,
		0,
		40.19,
		-3.00,
		-23.32,
		-13.11,
		-4.67,
		20194168.22,
		36.9892,
		105,
		0,
		40.19,
		-3.00,
		-23.32,
		26.31,
		-39.04,
		20671871.64,
		24.3477,
		105,
		4,
		40.19,
		-3.00,
		-23.32,
		79.33,
		-55.34,
		20679595.44,
		72.4723,
		105,
		4,
		40.19,
		-3.00,
		-23.32,
		107.19,
		-10.65,
		19943686.06,
		108.8940,
		105,
		4,
		40.19,
		-3.00,
		-23.32,
		56.35,
		47.54,
		20322471.38,
		66.2816,
		105,
		8,
		40.19,
		-3.00,
		-23.32,
		7.14,
		-19.55,
		20226657.34,
		98.5299,
		105,
		8,
		40.19,
		-3.00,
		-23.32,
		51.96,
		-1.90,
		20218595.37,
		78.9012,
		105,
		8,
		40.19,
		-3.00,
		-23.32,
		89.22,
		-40.56,
		20055109.63,
		149.9891,
		105,
		12,
		40.19,
		-3.00,
		-23.32,
		90.78,
		-28.26,
		20081398.25,
		153.7801,
		105,
		12,
		40.19,
		-3.00,
		-23.32,
		35.75,
		-14.88,
		20010521.91,
		88.5214,
		105,
		12,
		40.19,
		-3.00,
		-23.32,
		81.09,
		35.20,
		20278071.09,
		161.6302,
		105,
		16,
		40.19,
		-3.00,
		-23.32,
		14.89,
		32.88,
		20636367.52,
		86.2076,
		105,
		16,
		40.19,
		-3.00,
		-23.32,
		2.04,
		11.23,
		20394926.95,
		83.4856,
		105,
		16,
		40.19,
		-3.00,
		-23.32,
		22.79,
		-35.87,
		20125991.19,
		73.8076,
		105,
		20,
		40.19,
		-3.00,
		-23.32,
		54.11,
		3.15,
		20251696.28,
		36.3307,
		105,
		20,
		40.19,
		-3.00,
		-23.32,
		95.06,
		17.94,
		20246498.07,
		69.0422,
		105,
		20,
		40.19,
		-3.00,
		-23.32,
		-1.81,
		-52.00,
		20332764.38,
		69.9599,
		105,
		0,
		115.89,
		-31.80,
		12.78,
		119.90,
		-8.76,
		19941513.27,
		23.8358,
		105,
		0,
		115.89,
		-31.80,
		12.78,
		165.14,
		-13.93,
		20181976.57,
		38.5473,
		105,
		0,
		115.89,
		-31.80,
		12.78,
		76.65,
		-41.43,
		20157673.77,
		23.6849,
		105,
		4,
		115.89,
		-31.80,
		12.78,
		107.19,
		-10.65,
		19943685.24,
		40.8648,
		105,
		4,
		115.89,
		-31.80,
		12.78,
		79.33,
		-55.34,
		20679595.29,
		44.0008,
		105,
		4,
		115.89,
		-31.80,
		12.78,
		64.90,
		-17.58,
		20177185.06,
		57.6945,
		105,
		8,
		115.89,
		-31.80,
		12.78,
		127.35,
		23.46,
		19837695.71,
		82.6376,
		105,
		8,
		115.89,
		-31.80,
		12.78,
		89.22,
		-40.56,
		20055109.56,
		39.1752,
		105,
		8,
		115.89,
		-31.80,
		12.78,
		148.31,
		-29.93,
		20109263.99,
		40.6714,
		105,
		12,
		115.89,
		-31.80,
		12.78,
		90.78,
		-28.26,
		20081398.25,
		23.2595,
		105,
		12,
		115.89,
		-31.80,
		12.78,
		133.47,
		-24.87,
		19975574.41,
		21.5454,
		105,
		12,
		115.89,
		-31.80,
		12.78,
		166.97,
		-3.87,
		20196778.56,
		38.6379,
		105,
		16,
		115.89,
		-31.80,
		12.78,
		124.09,
		-14.31,
		20100697.90,
		15.0902,
		105,
		16,
		115.89,
		-31.80,
		12.78,
		154.31,
		-45.19,
		20116286.17,
		16.5133,
		105,
		16,
		115.89,
		-31.80,
		12.78,
		-167.50,
		-43.24,
		20095343.13,
		25.6441,
		105,
		20,
		115.89,
		-31.80,
		12.78,
		131.65,
		-31.56,
		20066111.12,
		7.5508,
		105,
		20,
		115.89,
		-31.80,
		12.78,
		115.68,
		-52.78,
		20231909.06,
		7.8470,
		105,
		20,
		115.89,
		-31.80,
		12.78,
		50.87,
		-50.69,
		20186511.77,
		12.4908,
	};
	
	double low[]=
	{
		105,
		0,
		141.13,
		39.14,
		117.00,
		165.14,
		-13.93,
		20181976.50,
		51.5270,
		105,
		0,
		141.13,
		39.14,
		117.00,
		85.59,
		36.64,
		20015444.79,
		25.3869,
		105,
		0,
		141.13,
		39.14,
		117.00,
		119.90,
		-8.76,
		19941513.27,
		40.6363,
		105,
		4,
		141.13,
		39.14,
		117.00,
		107.19,
		-10.65,
		19943685.88,
		67.1614,
		105,
		4,
		141.13,
		39.14,
		117.00,
		38.39,
		51.98,
		20457198.52,
		41.9432,
		105,
		4,
		141.13,
		39.14,
		117.00,
		-130.77,
		54.40,
		20121312.41,
		38.0286,
		105,
		8,
		141.13,
		39.14,
		117.00,
		179.50,
		51.35,
		19967933.94,
		20.6012,
		105,
		8,
		141.13,
		39.14,
		117.00,
		97.28,
		21.46,
		19941941.52,
		30.6540,
		105,
		8,
		141.13,
		39.14,
		117.00,
		84.26,
		54.68,
		20305726.98,
		25.3269,
		105,
		12,
		141.13,
		39.14,
		117.00,
		62.65,
		54.77,
		20370905.24,
		18.7723,
		105,
		12,
		141.13,
		39.14,
		117.00,
		115.63,
		-1.28,
		20165065.92,
		22.6981,
		105,
		12,
		141.13,
		39.14,
		117.00,
		81.09,
		35.20,
		20278071.22,
		19.0137,
		105,
		16,
		141.13,
		39.14,
		117.00,
		124.09,
		-14.31,
		20100698.19,
		20.1095,
		105,
		16,
		141.13,
		39.14,
		117.00,
		-130.60,
		49.21,
		20288020.98,
		12.3098,
		105,
		16,
		141.13,
		39.14,
		117.00,
		161.97,
		13.35,
		20265041.53,
		10.1019,
		105,
		20,
		141.13,
		39.14,
		117.00,
		84.18,
		36.59,
		19953853.27,
		10.2320,
		105,
		20,
		141.13,
		39.14,
		117.00,
		54.67,
		51.65,
		20511861.27,
		11.6913,
		105,
		20,
		141.13,
		39.14,
		117.00,
		-136.92,
		46.53,
		20309713.36,
		15.2671,
		105,
		0,
		-155.46,
		19.80,
		3754.69,
		165.14,
		-13.93,
		20181976.58,
		67.7247,
		105,
		0,
		-155.46,
		19.80,
		3754.69,
		179.32,
		9.92,
		20274303.54,
		41.9874,
		105,
		0,
		-155.46,
		19.80,
		3754.69,
		-144.16,
		-15.44,
		20007317.84,
		50.6599,
		105,
		4,
		-155.46,
		19.80,
		3754.69,
		-130.77,
		54.40,
		20121312.45,
		32.8732,
		105,
		4,
		-155.46,
		19.80,
		3754.69,
		-99.26,
		37.44,
		20066769.88,
		39.5971,
		105,
		4,
		-155.46,
		19.80,
		3754.69,
		-85.72,
		53.69,
		20544786.60,
		41.9149,
		105,
		8,
		-155.46,
		19.80,
		3754.69,
		178.35,
		-7.05,
		20372509.81,
		22.7673,
		105,
		8,
		-155.46,
		19.80,
		3754.69,
		-125.97,
		2.30,
		20251559.90,
		20.2296,
		105,
		8,
		-155.46,
		19.80,
		3754.69,
		179.50,
		51.35,
		19967934.29,
		19.9825,
		105,
		12,
		-155.46,
		19.80,
		3754.69,
		158.88,
		-12.61,
		20145417.20,
		22.1744,
		105,
		12,
		-155.46,
		19.80,
		3754.69,
		-146.53,
		22.03,
		20069033.97,
		9.7271,
		105,
		12,
		-155.46,
		19.80,
		3754.69,
		-153.30,
		-39.75,
		20672066.87,
		21.9415,
		105,
		16,
		-155.46,
		19.80,
		3754.69,
		-140.58,
		51.70,
		20455387.61,
		12.7349,
		105,
		16,
		-155.46,
		19.80,
		3754.69,
		-167.50,
		-43.24,
		20095343.11,
		25.0219,
		105,
		16,
		-155.46,
		19.80,
		3754.69,
		-164.50,
		27.08,
		20494802.61,
		10.4564,
		105,
		20,
		-155.46,
		19.80,
		3754.69,
		-172.71,
		-20.37,
		20225145.06,
		44.0358,
		105,
		20,
		-155.46,
		19.80,
		3754.69,
		-136.92,
		46.53,
		20309713.37,
		31.7940,
		105,
		20,
		-155.46,
		19.80,
		3754.69,
		-82.52,
		20.64,
		19937791.48,
		67.4750,
	};
	
	double error=0;

	double a0 = 236.831641;
	double a1 = -0.39362878;
	double a2 = 0.00402826613;
	
	for (int i=0 ; i<36 ; ++i)
	{
		int idx=i*9;
		
		double doy = high[idx + 0];
		double utc_time = high[idx + 1];
		double rxlon = pi * high[idx + 2] / 180.;
		double rxlat = pi * high[idx + 3] / 180.;
		double rxalt = high[idx + 4];
		double satlon = pi * high[idx + 5] / 180.;
		double satlat = pi * high[idx + 6] / 180.;
		double satalt = high[idx + 7];
		double output_stec = high[idx + 8];
		
		double STEC= stec(a0, a1, a2, rxlat, rxlon, rxalt, satlat, satlon, satalt, utc_time, doy);
		double diff = STEC - output_stec;
		error += diff * diff;
	}
	
	a0=121.129893;
	a1=0.351254133;
    a2=0.0134635348;
	
	for (int i=0 ; i<36 ; ++i)
	{
		int idx=i*9;
		
		double doy = med[idx + 0];
		double utc_time = med[idx + 1];
		double rxlon = pi * med[idx + 2] / 180.;
		double rxlat = pi * med[idx + 3] / 180.;
		double rxalt = med[idx + 4];
		double satlon = pi * med[idx + 5] / 180.;
		double satlat = pi * med[idx + 6] / 180.;
		double satalt = med[idx + 7];
		double output_stec = med[idx + 8];

		double STEC = stec(a0, a1, a2, rxlat, rxlon, rxalt, satlat, satlon, satalt, utc_time, doy);
		double diff = STEC - output_stec;
		error += diff * diff;
	}
	
	a0=2.580271;
	a1=0.127628236;
    a2=0.0252748384;
	
	for (int i=0 ; i<36 ; ++i)
	{
		int idx=i*9;
		
		double doy = low[idx + 0];
		double utc_time = low[idx + 1];
		double rxlon = pi * low[idx + 2] / 180.;
		double rxlat = pi * low[idx + 3] / 180.;
		double rxalt = low[idx + 4];
		double satlon = pi * low[idx + 5] / 180.;
		double satlat = pi * low[idx + 6] / 180.;
		double satalt = low[idx + 7];
		double output_stec = low[idx + 8];

		double STEC = stec(a0, a1, a2, rxlat, rxlon, rxalt, satlat, satlon, satalt, utc_time, doy);
		double diff = STEC - output_stec;
		error += diff*diff;
	}

	error = sqrt(error);
	
	return error;
}

};// namespace
