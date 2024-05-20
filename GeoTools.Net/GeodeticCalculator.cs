///////////////////////////////////////////////////////////////////////////////////
// adjusted code subset, ported from geotools version 2.3.3a.
//  -  Original source file: gt2-2.3.3\module\referencing\src\org\geotools\referencing\GeodeticCalculator.java
//
///////////////////////////////////////////////////////////////////////////////////
// Modifications of 1/30/2008:
// Retains all functionality of original Java/C# port.
//  - ComputeDirection() has been optimized for run-time.
//  - ArthmeticException throws have been replaced with Haversine distance call.
//  - New public method "Distance()" calls ComputeDirection directly.
// Modifications of 7/31/2008:
//  - Corrected problem using old methods:
//       SetStartingGeographicPoint(bosLat, bosLon);
//       SetDestinationGeographicPoint(laxLat, laxLon);
//       totalDistanceMeters = utilsGCalc.GetOrthodromicDistance();
///////////////////////////////////////////////////////////////////////////////////

/*
 *  Original license:
 *
 *    GeoTools - OpenSource mapping toolkit
 *    http://geotools.org
 *    (C) 2004-2006, GeoTools Project Managment Committee (PMC)
 *
 *    This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation;
 *    version 2.1 of the License.
 *
 *    This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    Portions of this file is adapted from Fortran code provided by NOAA.
 *    Programmed for CDC-6600 by LCDR L.Pfeifer NGS ROCKVILLE MD 18FEB75
 *    Modified for IBM SYSTEM 360 by John G.Gergen NGS ROCKVILLE MD 7507
 *    Source: ftp://ftp.ngs.noaa.gov/pub/pcsoft/for_inv.3d/source/
 */

namespace GeoTools.Net
{
    public sealed class GeodeticCalculator : IGeodeticCalculator
    {
        public GeodeticCalculator()
        {

            /* Calculation of GNPARC parameters */
            double E2 = eccentricitySquared;
            double E4 = E2 * E2;
            double E6 = E4 * E2;
            double E8 = E6 * E2;
            double EX = E8 * E2;

            A = 1.0 + 0.75 * E2 + 0.703125 * E4 + 0.68359375 * E6 + 0.67291259765625 * E8 + 0.6661834716796875 * EX;
            B = 0.75 * E2 + 0.9375 * E4 + 1.025390625 * E6 + 1.07666015625 * E8 + 1.1103057861328125 * EX;
            C = 0.234375 * E4 + 0.41015625 * E6 + 0.538330078125 * E8 + 0.63446044921875 * EX;
            D = 0.068359375 * E6 + 0.15380859375 * E8 + 0.23792266845703125 * EX;
            E = 0.01922607421875 * E8 + 0.0528717041015625 * EX;
            F = 0.00528717041015625 * EX;

            maxOrthodromicDistance = semiMajorAxis * (1.0 - E2) * Math.PI * A - 1.0;

            T1 = 1.0;
            T2 = -0.25 * f * (1.0 + f + f2);
            T4 = 0.1875 * f2 * (1.0 + 2.25 * f);
            T6 = 0.1953125 * f3;

            double a = f3 * (1.0 + 2.25 * f);
            a01 = -f2 * (1.0 + f + f2) / 4.0;
            a02 = 0.1875 * a;
            a03 = -0.1953125 * f4;
            a21 = -a01;
            a22 = -0.25 * a;
            a23 = 0.29296875 * f4;
            a42 = 0.03125 * a;
            a43 = 0.05859375 * f4;
            a63 = 5.0 * f4 / 768.0;

            /* calculation of GPNHRI parameters */
            f = (semiMajorAxis - semiMinorAxis) / semiMajorAxis;
            //f = (Constants.Constants.SEMI_MAJOR_AXIS - Constants.Constants.SEMI_MINOR_AXIS) / Constants.Constants.SEMI_MAJOR_AXIS;
            fo = 1.0 - f;
            f2 = f * f;
            f3 = f * f2;
            f4 = f * f3;
            // The eccentricity squared of the referenced ellipsoid.
            eccentricitySquared = f * (2.0 - f);
            ESQP = eccentricitySquared / (1.0 - eccentricitySquared);
            alimit = Math.PI * fo;

        }

        /**
           * Tolerance factors from the strictest (<CODE>TOLERANCE_0 </CODE>)
           * to the most relax one (<CODE>TOLERANCE_3</CODE>).
           */
        private const double TOLERANCE_0 = 5.0e-15, // tol0
            TOLERANCE_1 = 5.0e-14, // tol1
            TOLERANCE_2 = 5.0e-13, // tt
            TOLERANCE_3 = 7.0e-3; // tol2

        /*
         * The semi major axis of the referenced ellipsoid in meter
         */
        private const double semiMajorAxis = 6378137.0;

        /*
         * The semi minor axis of the referenced ellipsoid. WGS84
         */
        private const double semiMinorAxis = 6356752.314245179;

        private double f;
        private double fo;
        private double f2;
        private double f3;
        private double f4;

        // The eccentricity squared of the referenced ellipsoid.
        private double eccentricitySquared;
        private double ESQP;
        private double alimit;

        /*
         * The maximum orthodromic distance that could be calculated onto the referenced ellipsoid.
         */
        private double maxOrthodromicDistance;

        /**
         * GPNARC parameters computed from the ellipsoid.
         */
        private double A, B, C, D, E, F;

        /**
         * Parameters computed from the ellipsoid.
         */
        private double T1, T2, T4, T6;

        /**
         * Parameters computed from the ellipsoid.
         */
        private double a01, a02, a03, a21, a22, a23, a42, a43, a63;

        /**
         * The (<var>latitude</var>, <var>longitude</var>) coordinate of the first point
         * <strong>in radians</strong>. This point is set by {@link #setStartingGeographicPoint}.
         */
        private double lat1, long1;

        /**
         * The (<var>latitude</var>, <var>longitude</var>) coordinate of the destination point
         * <strong>in radians</strong>. This point is set by {@link #setDestinationGeographicPoint}.
         */
        private double lat2, long2;

        /**
         * The distance and azimuth (in radians) from the starting point
         * ({@link #long1}, {@link #lat1}) to the destination point
         * ({@link #long2}, {@link #lat2}).
         */
        private double distance, azimuth;

        /**
         * Tell if the destination point is valid.
         * {@code false} if {@link #long2} and {@link #lat2} need to be computed.
         */
        private bool destinationValid;

        /**
         * Tell if the azimuth and the distance are valid.
         * {@code false} if {@link #distance} and {@link #azimuth} need to be computed.
         */
        private bool directionValid;


        /// <summary>
        /// Set the starting point in geographic coordinates.
        /// The {@linkplain #getAzimuth() azimuth},
        /// the {@linkplain #getOrthodromicDistance() orthodromic distance} and
        /// the {@linkplain #getDestinationPoint() destination point} are discarded.
        /// They will need to be specified again.
        /// -- since 2.3
        /// </summary>
        /// <param name="latitude">The latitude  in decimal degrees between  -90&#x00B0; and  +90&#x00B0;</param>
        /// <param name="longitude">The longitude in decimal degrees between -180&#x00B0; and +180&#x00B0;</param>
        /// <throws>IllegalArgumentException if the longitude or the latitude is out of bounds.</throws>
        public void SetStartingGeographicPoint(double latitude, double longitude)
        {
            // Check first in case an exception is raised
            // (in other words, we change all or nothing).
            long1 = CheckLongitude(longitude);
            lat1 = CheckLatitude(latitude);

            destinationValid = false;
            directionValid = false;
        }

        /// <summary>
        /// Set the destination point in geographic coordinates. The azimuth and distance values
        /// will be updated as a side effect of this call. They will be recomputed the next time
        /// {@link #getAzimuth()} or {@link #getOrthodromicDistance()} are invoked.
        /// -- since 2.3
        /// </summary>
        /// <param name="latitude">The latitude in decimal degrees between  -90&#x00B0; and  +90&#x00B0;</param>
        /// <param name="longitude">The longitude in decimal degrees between  -90&#x00B0; and  +90&#x00B0;</param>
        /// <throws>IllegalArgumentException if the longitude or the latitude is out of bounds.</throws>
        public void SetDestinationGeographicPoint(double latitude, double longitude)
        {
            // Check first in case an exception is raised
            // (in other words, we change all or nothing).
            long2 = CheckLongitude(longitude);
            lat2 = CheckLatitude(latitude);

            destinationValid = true;
            directionValid = false;
        }

        public double GetDestinationLat()
        {
            if (!destinationValid)
            {
                ComputeDestinationPoint();
            }

            return AngleConverter.RadiansToDegrees(lat2);
        }

        public double GetDestinationLon()
        {
            if (!destinationValid)
            {
                ComputeDestinationPoint();
            }

            return AngleConverter.RadiansToDegrees(long2);
        }

        /// <summary>
        /// Set the azimuth and the distance from the {@linkplain #getStartingGeographicPoint
        /// starting point}. The destination point will be updated as a side effect of this call.
        /// It will be recomputed the next time {@link #getDestinationPoint()} is invoked.
        /// -- see #getAzimuth
        /// -- see #getOrthodromicDistance
        /// </summary>
        /// <param name="azimuth">The azimuth in decimal degrees from -180&#x00B0; to 180&#x00B0;.</param>
        /// <param name="distance">distance The orthodromic distance in the same units as the 
        /// {@linkplain #getEllipsoid ellipsoid} axis.</param>
        public void SetDirection(double azimuth, double distance)
        {
            // Check first in case an exception is raised
            // (in other words, we change all or nothing).
            azimuth = CheckAzimuth(azimuth); // degree => radian conversion
            //// checkOrthodromicDistance(distance);  cant do for now
            // Check passed. Now performs the changes in this object.
            this.azimuth = azimuth;
            this.distance = distance;
            destinationValid = false;
            directionValid = true;
        }

        /// <summary>
        /// Returns the azimuth. This method returns the value set by the last call to
        /// {@linkplain #setDirection(double,double) setDirection}(azimuth,distance),
        /// except if {@linkplain #setDestinationPoint(double,double) setDestinationPoint}(...)
        /// has been invoked after. In this later case, the azimuth will be computed from the
        /// {@linkplain #getStartingGeographicPoint starting point} to the destination point.
        /// </summary>
        /// <returns>The azimuth, in decimal degrees from -180&#x00B0; to +180&#x00B0;.</returns>
        /// <throws>IllegalStateException if the destination point has not been set.</throws>
        public double GetAzimuth()
        {
            if (!directionValid)
            {
                ComputeDirection();
            }

            return AngleConverter.RadiansToDegrees(azimuth);
        }

        /// <summary>
        /// Returns the orthodromic distance. This method returns the value set by the last call to
        /// {@linkplain #setDirection(double,double) setDirection}(azimuth,distance),
        /// except if {@linkplain #setDestinationPoint(double,double) setDestinationPoint}(...)
        /// has been invoked after. In this later case, the distance will be computed from the
        /// {@linkplain #getStartingGeographicPoint starting point} to the destination point.
        /// </summary>
        /// <returns>The orthodromic distance, in the same units as the 
        /// {@linkplain #getEllipsoid ellipsoid} axis.</returns>
        /// <throws>IllegalStateException if the destination point has not been set.</throws>
        public double GetOrthodromicDistance()
        {
            if (!directionValid)
            {
                ComputeDirection();
            }

            return distance;
        }

        /// <summary>
        /// Returns distance in meters.
        /// </summary>
        /// <param name="latitude1">Degrees of latitude</param>
        /// <param name="longitude1">Degrees of longitude</param>
        /// <param name="latitude2">Degrees of latitude</param>
        /// <param name="longitude2">Degrees of longitude</param>
        /// <returns>distance in meters</returns>
        public double Distance(double latitude1, double longitude1, double latitude2, double longitude2)
        {
            long1 = AngleConverter.DegreesToRadians(longitude1);
            lat1 = AngleConverter.DegreesToRadians(latitude1);

            long2 = AngleConverter.DegreesToRadians(longitude2);
            lat2 = AngleConverter.DegreesToRadians(latitude2);

            ComputeDirection();

            return distance;
        }

        /// <summary>
        /// Optimization of original HaversineFormulaGreatCircleDistanceInKm().
        /// </summary>
        /// <returns></returns>
        private double HaversineGCDistanceMeters()
        {

            /*
            The Haversine formula according to Dr. Math.
            http://mathforum.org/library/drmath/view/51879.html

            dlon = lon2 - lon1
            dlat = lat2 - lat1
            a = (sin(dlat/2))^2 + cos(lat1) * cos(lat2) * (sin(dlon/2))^2
            c = 2 * atan2(sqrt(a), sqrt(1-a))
            d = R * c

            Where
                * dlon is the change in longitude
                * dlat is the change in latitude
                * c is the great circle distance in Radians.
                * R is the radius of a spherical Earth.
                * The locations of the two points in
                    spherical coordinates (longitude and
                    latitude) are lon1,lat1 and lon2, lat2.
            */

            double dLat1InRad = lat1;
            double dLong1InRad = long1;
            double dLat2InRad = lat2;
            double dLong2InRad = long2;

            double dLongitude = dLong2InRad - dLong1InRad;
            double dLatitude = dLat2InRad - dLat1InRad;

            // Intermediate result a.
            double a = Math.Pow(Math.Sin(dLatitude / 2.0), 2.0) +
                       Math.Cos(dLat1InRad) * Math.Cos(dLat2InRad) *
                       Math.Pow(Math.Sin(dLongitude / 2.0), 2.0);

            // Intermediate result c (great circle distance in Radians).
            double c = 2.0 * Math.Atan2(Math.Sqrt(a), Math.Sqrt(1.0 - a));

            // Distance.
            const Double kEarthRadiusMeters = 6376500;
            double dDistance = kEarthRadiusMeters * c;

            return dDistance;
        }

        /**
       * Returns an angle between -{@linkplain Math#PI PI} and {@linkplain Math#PI PI}
       * equivalent to the specified angle in radians.
       *
       * @param  alpha An angle value in radians.
       * @return The angle between -{@linkplain Math#PI PI} and {@linkplain Math#PI PI}.
       *
       */
        private static double CastToAngleRange(double alpha)
        {
            return alpha - (2 * Math.PI) * Math.Floor(alpha / (2 * Math.PI) + 0.5);
        }

        /**
         * Checks the latitude validity. The argument {@code latitude} should be
         * greater or equal than -90 degrees and lower or equals than +90 degrees. As
         * a convenience, this method returns the latitude in radians.
         *
         * @param  latitude The latitude value in <strong>decimal degrees</strong>.
         * @return The latitude value in <strong>radians</strong>.
         * @throws IllegalArgumentException if {@code latitude} is not between -90 and +90 degrees.
         */
        private static double CheckLatitude(double latitude)
        {
            if (latitude >= -90.0 && latitude <= 90.0)
            {
                return AngleConverter.DegreesToRadians(latitude);
            }

            throw new ArgumentOutOfRangeException(nameof(latitude));
        }

        /**
         * Checks the longitude validity. The argument {@code longitude} should be
         * greater or equal than -180 degrees and lower or equals than +180 degrees. As
         * a convenience, this method returns the longitude in radians.
         *
         * @param  longitude The longitude value in <strong>decimal degrees</strong>.
         * @return The longitude value in <strong>radians</strong>.
         * @throws IllegalArgumentException if {@code longitude} is not between -180 and +180 degrees.
         */
        private static double CheckLongitude(double longitude)
        {
            if (longitude >= -180.0 && longitude <= 180.0)
            {
                return AngleConverter.DegreesToRadians(longitude);
            }

            throw new ArgumentOutOfRangeException(nameof(longitude));
        }

        /**
         * Checks the azimuth validity. The argument {@code azimuth}  should be
         * greater or equal than -180 degrees and lower or equals than +180 degrees.
         * As a convenience, this method returns the azimuth in radians.
         *
         * @param  azimuth The azimuth value in <strong>decimal degrees</strong>.
         * @return The azimuth value in <strong>radians</strong>.
         * @throws IllegalArgumentException if {@code azimuth} is not between -180 and +180 degrees.
         */
        private static double CheckAzimuth(double azimuth)
        {
            if (azimuth >= -180.0 && azimuth <= 180.0)
            {
                return AngleConverter.DegreesToRadians(azimuth);
            }

            throw new ArgumentOutOfRangeException(nameof(azimuth));
        }

        /**
         * Checks the orthodromic distance validity. Arguments {@code orthodromicDistance}
         * should be greater or equal than 0 and lower or equals than the maximum orthodromic distance.
         *
         * @param  distance The orthodromic distance value.
         * @throws IllegalArgumentException if {@code orthodromic distance} is not between
         *                                  0 and the maximum orthodromic distance.
         */
        private void CheckOrthodromicDistance(double distance)
        {
            if (!(distance >= 0.0 && distance <= maxOrthodromicDistance))
            {
                throw new ArgumentOutOfRangeException(nameof(distance));
            }
        }

        /**
         * Computes the destination point from the {@linkplain #getStartingGeographicPoint starting
         * point}, the {@linkplain #getAzimuth azimuth} and the {@linkplain #getOrthodromicDistance
         * orthodromic distance}.
         *
         * @throws IllegalStateException if the azimuth and the distance have not been set.
         *
         * @see #getDestinationPoint
         */
        private void ComputeDestinationPoint()
        {
            if (!directionValid)
            {
                throw new Exception();
            }

            // Protect internal variables from changes
            double lat1 = this.lat1;
            double long1 = this.long1;
            double azimuth = this.azimuth;
            double distance = this.distance;

            /*
             * Solution of the geodetic direct problem after T.Vincenty.
             * Modified Rainsford's method with Helmert's elliptical terms.
             * Effective in any azimuth and at any distance short of antipodal.
             *
             * Latitudes and longitudes in radians positive North and East.
             * Forward azimuths at both points returned in radians from North.
             *
             * Programmed for CDC-6600 by LCDR L.Pfeifer NGS ROCKVILLE MD 18FEB75
             * Modified for IBM SYSTEM 360 by John G.Gergen NGS ROCKVILLE MD 7507
             * Ported from Fortran to Java by Daniele Franzoni.
             *
             * Source:
             * *         subroutine DIRECT1
             */
            double TU = fo * Math.Sin(lat1) / Math.Cos(lat1);
            double SF = Math.Sin(azimuth);
            double CF = Math.Cos(azimuth);
            double BAZ = (CF != 0) ? Math.Atan2(TU, CF) * 2.0 : 0;
            double CU = 1 / Math.Sqrt(TU * TU + 1.0);
            double SU = TU * CU;
            double SA = CU * SF;
            double C2A = 1.0 - SA * SA;
            double X = Math.Sqrt((1.0 / fo / fo - 1) * C2A + 1.0) + 1.0;
            X = (X - 2.0) / X;
            double C = 1.0 - X;
            C = (X * X / 4.0 + 1.0) / C;
            double D = (0.375 * X * X - 1.0) * X;
            TU = distance / fo / semiMajorAxis / C;
            double Y = TU;
            double SY, CY, CZ, E;
            do
            {
                SY = Math.Sin(Y);
                CY = Math.Cos(Y);
                CZ = Math.Cos(BAZ + Y);
                E = CZ * CZ * 2.0 - 1.0;
                C = Y;
                X = E * CY;
                Y = E + E - 1.0;
                Y = (((SY * SY * 4.0 - 3.0) * Y * CZ * D / 6.0 + X) * D / 4.0 - CZ) * SY * D + TU;
            } while (Math.Abs(Y - C) > TOLERANCE_1);

            BAZ = CU * CY * CF - SU * SY;
            C = fo * Math.Sqrt(SA * SA + BAZ * BAZ);
            D = SU * CY + CU * SY * CF;
            lat2 = Math.Atan2(D, C);
            C = CU * CY - SU * SY * CF;
            X = Math.Atan2(SY * SF, C);
            C = ((-3.0 * C2A + 4.0) * f + 4.0) * C2A * f / 16.0;
            D = ((E * CY * C + CZ) * SY * C + Y) * SA;
            long2 = long1 + X - (1.0 - C) * D * f;
            long2 = CastToAngleRange(long2);
            destinationValid = true;
        }

        /**
         * Calculates the meridian arc length between two points in the same meridian
         * in the referenced ellipsoid.
         *
         * @param  P1 The latitude of the first  point (in radians).
         * @param  P2 The latitude of the second point (in radians).
         * @return Returned the meridian arc length between P1 and P2
         */
        private double GetMeridianArcLengthRadians(double P1, double P2)
        {
            /*
             * Latitudes P1 and P2 in radians positive North and East.
             * Forward azimuths at both points returned in radians from North.
             *
             * Source: ftp://ftp.ngs.noaa.gov/pub/pcsoft/for_inv.3d/source/inverse.for
             *         subroutine GPNARC
             *         version    200005.26
             *         written by Robert (Sid) Safford
             *
             * Ported from Fortran to Java by Daniele Franzoni.
             */
            double S1 = Math.Abs(P1);
            double S2 = Math.Abs(P2);
            double DA = (P2 - P1);

            // Check for a 90 degree lookup
            if (S1 > TOLERANCE_0 || S2 <= (Math.PI / 2 - TOLERANCE_0) || S2 >= (Math.PI / 2 + TOLERANCE_0))
            {
                double DB = Math.Sin(P2 * 2.0) - Math.Sin(P1 * 2.0);
                double DC = Math.Sin(P2 * 4.0) - Math.Sin(P1 * 4.0);
                double DD = Math.Sin(P2 * 6.0) - Math.Sin(P1 * 6.0);
                double DE = Math.Sin(P2 * 8.0) - Math.Sin(P1 * 8.0);
                double DF = Math.Sin(P2 * 10.0) - Math.Sin(P1 * 10.0);

                // Compute the S2 part of the series expansion
                S2 = -DB * B / 2.0 + DC * C / 4.0 - DD * D / 6.0 + DE * E / 8.0 - DF * F / 10.0;
            }

            // Compute the S1 part of the series expansion
            S1 = DA * A;

            // Compute the arc length
            return Math.Abs(semiMajorAxis * (1.0 - eccentricitySquared) * (S1 + S2));
        }

        /**
         * Computes the azimuth and orthodromic distance from the
         * {@linkplain #getStartingGeographicPoint starting point} and the
         * {@linkplain #getDestinationPoint destination point}.
         *
         * @throws IllegalStateException if the destination point has not been set.
         *
         * @see #getAzimuth
         * @see #getOrthodromicDistance
         */
        private void ComputeDirection()
        {
            /*
             * Solution of the geodetic inverse problem after T.Vincenty.
             * Modified Rainsford's method with Helmert's elliptical terms.
             * Effective in any azimuth and at any distance short of antipodal.
             *
             * Latitudes and longitudes in radians positive North and East.
             * Forward azimuths at both points returned in radians from North.
             *
             * Programmed for CDC-6600 by LCDR L.Pfeifer NGS ROCKVILLE MD 18FEB75
             * Modified for IBM SYSTEM 360 by John G.Gergen NGS ROCKVILLE MD 7507
             * Ported from Fortran to Java by Daniele Franzoni.
             *
             * Source: ftp://ftp.ngs.noaa.gov/pub/pcsoft/for_inv.3d/source/inverse.for
             *         subroutine GPNHRI
             *         version    200208.09
             *         written by robert (sid) safford
             */
            double dlon = CastToAngleRange(long2 - long1);
            double ss = Math.Abs(dlon);

            if (ss < TOLERANCE_1)
            {
                distance = GetMeridianArcLengthRadians(lat1, lat2);
                azimuth = (lat2 > lat1) ? 0.0 : Math.PI;
                directionValid = true;
                return;
            }

            /*
             * Computes the limit in longitude (alimit), it is equal
             * to twice  the distance from the equator to the pole,
             * as measured along the equator
             */

            if (ss >= alimit &&
                lat1 < TOLERANCE_3 && lat1 > -TOLERANCE_3 &&
                lat2 < TOLERANCE_3 && lat2 > -TOLERANCE_3)
            {
                // Computes an approximate AZ
                double CONS = (Math.PI - ss) / (Math.PI * f);
                double AZ = Math.Asin(CONS);
                double AZ_TEMP, S, AO;
                int iter = 0;
                do
                {
                    if (++iter > 8)
                    {
                        distance = HaversineGCDistanceMeters();
                        return;
                    }

                    S = Math.Cos(AZ);
                    double C2 = S * S;
                    // Compute new AO
                    AO = T1 + T2 * C2 + T4 * C2 * C2 + T6 * C2 * C2 * C2;
                    double CS = CONS / AO;
                    S = Math.Asin(CS);
                    AZ_TEMP = AZ;
                    AZ = S;
                } while (Math.Abs(S - AZ_TEMP) >= TOLERANCE_2);

                double AZ1 = (dlon < 0.0) ? 2.0 * Math.PI - S : S;
                azimuth = CastToAngleRange(AZ1);
                double AZ2 = 2.0 * Math.PI - AZ1;
                S = Math.Cos(AZ1);

                // Equatorial - geodesic(S-s) SMS
                double U2 = ESQP * S * S;
                double U4 = U2 * U2;
                double U6 = U4 * U2;
                double U8 = U6 * U2;
                double BO = 1.0 +
                            0.25 * U2 +
                            0.046875 * U4 +
                            0.01953125 * U6 +
                            -0.01068115234375 * U8;
                S = Math.Sin(AZ1);
                double SMS = semiMajorAxis * Math.PI * (1.0 - f * Math.Abs(S) * AO - BO * fo);
                distance = semiMajorAxis * ss - SMS;
                directionValid = true;
                return;
            }

            // the reduced latitudes
            double u1 = Math.Atan(fo * Math.Sin(lat1) / Math.Cos(lat1));
            double u2 = Math.Atan(fo * Math.Sin(lat2) / Math.Cos(lat2));
            double su1 = Math.Sin(u1);
            double cu1 = Math.Cos(u1);
            double su2 = Math.Sin(u2);
            double cu2 = Math.Cos(u2);
            double xy, w, q2, q4, q6, r2, r3, sig, ssig, slon, clon, sinalf, ab = dlon;
            int kcount = 0;

            do
            {
                if (++kcount > 8)
                {
                    distance = HaversineGCDistanceMeters();
                    return;
                }

                clon = Math.Cos(ab);
                slon = Math.Sin(ab);
                double csig = su1 * su2 + cu1 * cu2 * clon;
                ssig = Math.Sqrt(slon * cu2 * slon * cu2 +
                                 (su2 * cu1 - su1 * cu2 * clon) * (su2 * cu1 - su1 * cu2 * clon));

                sig = Math.Atan2(ssig, csig);
                sinalf = cu1 * cu2 * slon / ssig;
                w = (1.0 - sinalf * sinalf);
                double t4 = w * w;
                double t6 = w * t4;

                // the coefficients of type a
                double ao = f + a01 * w + a02 * t4 + a03 * t6;

                double a2 = a21 * w + a22 * t4 + a23 * t6;
                double a4 = a42 * t4 + a43 * t6;
                double a6 = a63 * t6;

                // the multiple angle functions
                double qo = 0.0;
                if (w > TOLERANCE_0)
                {
                    qo = -2.0 * su1 * su2 / w;
                    //qo = -2.0 * su1Xsu2 / w;
                }

                q2 = csig + qo;
                q4 = 2.0 * q2 * q2 - 1.0;
                q6 = q2 * (4.0 * q2 * q2 - 3.0);
                r2 = 2.0 * ssig * csig;
                r3 = ssig * (3.0 - 4.0 * ssig * ssig);

                // the longitude difference
                double s = sinalf * (ao * sig + a2 * ssig * q2 + a4 * r2 * q4 + a6 * r3 * q6);
                double xz = dlon + s;
                xy = Math.Abs(xz - ab);
                ab = dlon + s;
            } while (xy >= TOLERANCE_1);

            double z = ESQP * w;
            double bo = 1.0 + z * (1.0 / 4.0 + z * (-3.0 / 64.0 + z * (5.0 / 256.0 - z * (175.0 / 16384.0))));
            double b2 = z * (-1.0 / 4.0 + z * (1.0 / 16.0 + z * (-15.0 / 512.0 + z * (35.0 / 2048.0))));
            double b4 = z * z * (-1.0 / 128.0 + z * (3.0 / 512.0 - z * (35.0 / 8192.0)));
            double b6 = z * z * z * (-1.0 / 1536.0 + z * (5.0 / 6144.0));

            // The distance in ellipsoid axis units.
            distance = semiMinorAxis * (bo * sig + b2 * ssig * q2 + b4 * r2 * q4 + b6 * r3 * q6);
            double az1 = (dlon < 0.0) ? Math.PI * (3.0 / 2.0) : Math.PI / 2.0;

            // now compute the az1 & az2 for latitudes not on the equator
            if ((Math.Abs(su1) >= TOLERANCE_0) || (Math.Abs(su2) >= TOLERANCE_0))
            {
                double tana1 = slon * cu2 / (su2 * cu1 - clon * su1 * cu2);
                double sina1 = sinalf / cu1;

                // azimuths from north,longitudes positive east  
                az1 = Math.Atan2(sina1, sina1 / tana1);
            }

            azimuth = CastToAngleRange(az1);
            directionValid = true;
        }
    }
}