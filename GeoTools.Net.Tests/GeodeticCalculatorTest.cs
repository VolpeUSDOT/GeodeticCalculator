using FluentAssertions;

namespace GeoTools.Net.Tests
{
    [TestClass]
    public class GeodeticCalculatorTest
    {
        [Description("Tests GeodeticCalculator.Distance function")]
        [TestMethod]
        public void TestDistance()
        {
            // start and end point set to calculate geodetic distance
            GeodeticCalculator calc = new GeodeticCalculator();
            double bosLaDist = calc.Distance(42.36, -71.00, 33.94, -118.40);
            bosLaDist.Should().BeGreaterThan(4200000);
        }

        [Description("Tests Destination in LatLon over Distance from GeodeticCalculator")]
        [TestMethod]
        public void TestGetDestinationLatLon()
        {
            // start and end point set to calculate geodetic distance
            GeodeticCalculator calc = new GeodeticCalculator();
            calc.Distance(42.36, -71.00, 33.94, -118.40);
            double lat = calc.GetDestinationLat();
            double lon = calc.GetDestinationLon();
            lat.Should().BeApproximately(33.94, 0.0001);
            lon.Should().BeApproximately(-118.40, 0.0001);
        }

        [Description("Tests 0.0001 precision in LatLon over getting Destination from GeodeticCalculator")]
        [TestMethod]
        public void TestSetGeographicPoints()
        {
            GeodeticCalculator calc = new GeodeticCalculator();
            calc.SetStartingGeographicPoint(42.36, -71.00);
            calc.SetDestinationGeographicPoint(33.94, -118.40);
            double lat = calc.GetDestinationLat();
            double lon = calc.GetDestinationLon();

            lat.Should().BeApproximately(33.94, 0.0001);
            lon.Should().BeApproximately(-118.40, 0.0001);
        }

        [Description("Tests GetAzimuth from GeodeticCalculator with and without tolerance")]
        [TestMethod]
        public void TestGetAzimuth()
        {
            // set start and dest to calculate azimuth
            GeodeticCalculator calc = new GeodeticCalculator();
            calc.SetStartingGeographicPoint(42.0, -71.0);

            calc.SetDestinationGeographicPoint(43.0, -71.0);
            double azimuth = calc.GetAzimuth();
            azimuth.Should().Be(0.0);

            calc.SetDestinationGeographicPoint(41.0, -71.0);
            azimuth = calc.GetAzimuth();
            azimuth.Should().Be(180.0);

            double tolerance = 4.0; // IsWithInTolerance value

            calc.SetDestinationGeographicPoint(42.0, -72.0);
            azimuth = calc.GetAzimuth(); // -89.66542952260451
            azimuth.Should().BeApproximately(-90, tolerance);

            calc.SetDestinationGeographicPoint(42.33, -70.5);
            azimuth = calc.GetAzimuth(); // 48.254754022065534,
            azimuth.Should().BeApproximately(45.0, tolerance);

            calc.SetDestinationGeographicPoint(42.0, -70.5);
            azimuth = calc.GetAzimuth(); // 89.832716526744065, 
            azimuth.Should().BeApproximately(90.0, tolerance);

            calc.SetDestinationGeographicPoint(41.67, -70.5);
            azimuth = calc.GetAzimuth(); // 131.26212366423516
            azimuth.Should().BeApproximately(135.0, tolerance);

            calc.SetDestinationGeographicPoint(41.67, -71.5);
            azimuth = calc.GetAzimuth(); // -131.2621366423445
            azimuth.Should().BeApproximately(-135.0, tolerance);

            calc.SetDestinationGeographicPoint(42.33, -71.5);
            azimuth = calc.GetAzimuth(); // -48.25474022066258
            azimuth.Should().BeApproximately(-45.0, tolerance);
        }

        [Description("Tests Destination using Azimuth from GeodeticCalculator")]
        [TestMethod]
        public void TestDestinationUsingAzimuth()
        {
            TestDestination(42.0, -71.0, 43.0, -71.0);
            TestDestination(42.0, -71.0, 41.0, -71.0);
            TestDestination(42.0, -71.0, 42.0, -72.0);
            TestDestination(42.0, -71.0, 42.33, -70.5);
            TestDestination(42.0, -71.0, 42.0, -70.5);
            TestDestination(42.0, -71.0, 41.67, -70.5);
            TestDestination(42.0, -71.0, 41.67, -71.5);
            TestDestination(42.0, -71.0, 42.33, -71.5);
        }

        private void TestDestination(double startLat, double startLon, double destLat, double destLon)
        {
            GeodeticCalculator geoCalc = new GeodeticCalculator();
            geoCalc.SetStartingGeographicPoint(startLat, startLon);
            geoCalc.SetDestinationGeographicPoint(destLat, destLon);
            double azimuth = geoCalc.GetAzimuth();
            double distance = geoCalc.GetOrthodromicDistance();
            geoCalc.SetDirection(azimuth, distance);
            double destLatNew = geoCalc.GetDestinationLat();
            double destLonNew = geoCalc.GetDestinationLon();

            Console.WriteLine($"DiffLat: {destLatNew} vs. {destLat}, delta {destLatNew - destLat}");
            destLatNew.Should().BeApproximately(destLat, 0.005);
            Console.WriteLine($"DiffLon: {destLonNew} vs. {destLon}, delta {destLonNew - destLon}\n");
            destLonNew.Should().BeApproximately(destLon, 0.005);
        }
    }
}
