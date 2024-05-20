namespace GeoTools.Net
{
    public static class AngleConverter
    {
        /// <summary>
        /// Method to convert radians to degrees
        /// </summary>
        public static double RadiansToDegrees(double radians)
        {
            return radians * (180 / Math.PI);
        }

        /// <summary>
        /// Method to convert degrees to radians
        /// </summary>
        public static double DegreesToRadians(double degrees)
        {
            return degrees * (Math.PI / 180);
        }
    }
}
