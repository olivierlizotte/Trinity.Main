using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Trinity
{
    public class ElutionCurve
    {
        public double Area = 0.0;
        public double[] Coefficients = null;        

        public List<double> time = null;
        public List<double> intensityCount = null;
        public static ElutionCurve Create(PeptideSpectrumMatches psmMatches)
        {
            ElutionCurve theCurve = new ElutionCurve();
            // -- Test curve fitting function -- //
            theCurve.time = new List<double>();
            theCurve.intensityCount = new List<double>();
            foreach (PeptideSpectrumMatch psm in psmMatches)
            {
                theCurve.time.Add(psm.Query.spectrum.RetentionTimeInMin * 60 * 1000);
                theCurve.intensityCount.Add(psm.Query.spectrum.PrecursorIntensityPerMilliSecond);
            }
            theCurve.Compute();
            return theCurve;
        }

        public static ElutionCurve Create(List<Query> queries)
        {
            ElutionCurve theCurve = new ElutionCurve();
            // -- Test curve fitting function -- //
            theCurve.time = new List<double>();
            theCurve.intensityCount = new List<double>();
            foreach(Query query in queries)
            {
                theCurve.time.Add(query.spectrum.RetentionTimeInMin * 60 * 1000);
                theCurve.intensityCount.Add(query.spectrum.PrecursorIntensityPerMilliSecond);
            }
            theCurve.Compute();
            return theCurve;
        }

        public static ElutionCurve Create(List<double> intensitiesPerMilliSeconds, List<double> timesInMilliSeconds)
        {
            ElutionCurve theCurve = new ElutionCurve();
            // -- Test curve fitting function -- //
            theCurve.time = new List<double>();
            theCurve.intensityCount = new List<double>();
            for (int i = 0; i < intensitiesPerMilliSeconds.Count; i++)
            {
                theCurve.time.Add(timesInMilliSeconds[i]);
                theCurve.intensityCount.Add(intensitiesPerMilliSeconds[i]);
            }
            theCurve.Compute();
            return theCurve;
        }

        public static ElutionCurve Create(List<ProductSpectrum> spectras)
        {
            ElutionCurve theCurve = new ElutionCurve();
            // -- Test curve fitting function -- //
            theCurve.time = new List<double>();
            theCurve.intensityCount = new List<double>();
            foreach (ProductSpectrum spec in spectras)
            {
                theCurve.time.Add(spec.RetentionTimeInMin * 60 * 1000);
                theCurve.intensityCount.Add(spec.PrecursorIntensityPerMilliSecond);
            }
            theCurve.Compute();
            return theCurve;
        }

        public void Compute()
        {
            if (time != null && time.Count > 2)
            {
                double area1 = Proteomics.Utilities.CurveFitter.FitToPolynomial(time.ToArray(), intensityCount.ToArray(), out Coefficients);
                //double areaInterpol = Proteomics.Utilities.CurveFitter.AreaUnderTheCurve(time, intensityCount);
                double area2 = Proteomics.Utilities.CurveFitter.AreaUnderTheCurve(time, intensityCount);
                if(area1 / area2 > 2 || area1 / area2 < 0.5)
                    Console.WriteLine("Too much diff");
                Area = area1;// (area1 + area2) * 0.5;
                //if (Area > areaInterpol * 1.2 || Area < areaInterpol * 0.8)
                //    Console.WriteLine("Too much diff");
            }
            else
            {
                Area = 0;
                Coefficients = new double[0];
            }
        }

        public void AddPoint(double newTimePoint, double newIntensityPerMilliSeconds)
        {
            if (time == null || intensityCount == null)
            {
                time = new List<double>();
                intensityCount = new List<double>();
            }
            time.Add(newTimePoint);
            intensityCount.Add(newIntensityPerMilliSeconds);
        }
    }
}
