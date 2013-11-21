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

        public void Compute()
        {
            if (time != null && time.Count > 2)
                Area = Proteomics.Utilities.CurveFitter.FitToPolynomial(time.ToArray(), intensityCount.ToArray(), out Coefficients);
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
