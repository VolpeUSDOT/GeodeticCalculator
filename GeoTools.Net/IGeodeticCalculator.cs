namespace GeoTools.Net
{
    public interface IGeodeticCalculator
    {
        /// <summary>
        /// Set the starting point in geographic coordinates.
        /// The {@linkplain #getAzimuth() azimuth},
        /// the {@linkplain #getOrthodromicDistance() orthodromic distance} and
        /// the {@linkplain #getDestinationPoint() destination point} are discarded.
        /// </summary>
        /// <param name="latitude">The latitude  in decimal degrees between  -90&#x00B0; and  +90&#x00B0;</param>
        /// <param name="longitude">The longitude in decimal degrees between -180&#x00B0; and +180&#x00B0;</param>
        void SetStartingGeographicPoint(double latitude, double longitude);

        /// <summary>
        /// Set the destination point in geographic coordinates. The azimuth and distance values
        /// will be updated as a side effect of this call. They will be recomputed the next time
        /// {@link #getAzimuth()} or {@link #getOrthodromicDistance()} are invoked.
        /// </summary>
        /// <param name="latitude">The latitude in decimal degrees between  -90&#x00B0; and  +90&#x00B0;</param>
        /// <param name="longitude">The longitude in decimal degrees between  -90&#x00B0; and  +90&#x00B0;</param>
        void SetDestinationGeographicPoint(double latitude, double longitude);

        /// <summary>
        /// GetDestinationLat
        /// </summary>
        /// <returns>latitude in degrees</returns>
        double GetDestinationLat();

        /// <summary>
        /// GetDestinationLongitude
        /// </summary>
        /// <returns>longitude in degrees</returns>
        double GetDestinationLon();

        /// <summary>
        /// Set the azimuth and the distance from the {@linkplain #getStartingGeographicPoint
        /// starting point}. The destination point will be updated as a side effect of this call.
        /// It will be recomputed the next time {@link #getDestinationPoint()} is invoked.
        /// -- see #getAzimuth
        /// -- see #getOrthodromicDistance
        /// </summary>
        /// <param name="azimuth">The azimuth in decimal degrees from -180&#x00B0; to 180&#x00B0;.</param>
        /// <param name="distance">distance The orthodromic distance in the same units as the 
        /// {@linkplain #getEllipsoid ellipsoid} axis, i.e. in meters</param>
        void SetDirection(double azimuth, double distance);

        /// <summary>
        /// Returns the azimuth. This method returns the value set by the last call to
        /// {@linkplain #setDirection(double,double) setDirection}(azimuth,distance),
        /// except if {@linkplain #setDestinationPoint(double,double) setDestinationPoint}(...)
        /// has been invoked after. In this later case, the azimuth will be computed from the
        /// {@linkplain #getStartingGeographicPoint starting point} to the destination point.
        /// </summary>
        /// <returns>The azimuth, in decimal degrees from -180&#x00B0; to +180&#x00B0;.</returns>
        double GetAzimuth();

        /// <summary>
        /// Returns the orthodromic distance. This method returns the value set by the last call to
        /// {@linkplain #setDirection(double,double) setDirection}(azimuth,distance),
        /// except if {@linkplain #setDestinationPoint(double,double) setDestinationPoint}(...)
        /// has been invoked after. In this later case, the distance will be computed from the
        /// {@linkplain #getStartingGeographicPoint starting point} to the destination point.
        /// </summary>
        /// <returns>The orthodromic distance, in the same units as the 
        /// {@linkplain #getEllipsoid ellipsoid} axis, i.e. meters</returns>
        double GetOrthodromicDistance();

        /// <summary>
        /// Returns orthodromic distance in meters.
        /// </summary>
        /// <param name="latitude1">Degrees of latitude</param>
        /// <param name="longitude1">Degrees of longitude</param>
        /// <param name="latitude2">Degrees of latitude</param>
        /// <param name="longitude2">Degrees of longitude</param>
        /// <returns>distance in meters</returns>
        double Distance(double latitude1, double longitude1, double latitude2, double longitude2);
    }
}
