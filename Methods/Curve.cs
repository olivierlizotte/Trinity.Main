using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Trinity
{
    public class MaxFlowElutionCurve : ElutionCurve
    {
        public int nbProducts;
        public MaxFlowElutionCurve(int nbProductsUsed)
        {
            this.nbProducts = nbProductsUsed;
        }
    }

    public class ElutionCurveMerger
    {
        private List<ElutionCurve> Curves = new List<ElutionCurve>();
        private List<double>       Factor = new List<double>();
        public void AddCurve(ElutionCurve newCurve, double weight)
        {
            Curves.Add(newCurve);
            Factor.Add(weight);
        }

        public ElutionCurve Merge()
        {
            if(Curves.Count > 1)
            {
                ElutionCurve newCurve = new ElutionCurve();
                double sum = 0.0;
                foreach (double val in Factor)
                    sum += val;

                Dictionary<double, int> times = new Dictionary<double, int>();
                foreach (ElutionCurve curve in Curves)
                    foreach (double timePoint in curve.time)
                        if (!times.ContainsKey(timePoint))
                            times.Add(timePoint, 1);
                        else
                            times[timePoint]++;

                foreach(double key in times.Keys)
                    if (times[key] > 1)
                    {
                        double cumulIntensity = 0.0;
                        for (int i = 0; i < Curves.Count; i++)
                            cumulIntensity += Curves[i].InterpolateIntensity(key) * Factor[i] / sum;

                        newCurve.AddPoint(key, cumulIntensity);
                    }
                return newCurve;
            }
            else if (Curves.Count == 1)
                return Curves[0];
            return new ElutionCurve();
        }
    }

    public class ElutionCurve
    {
        public double Area = 0.0;
        public double[] Coefficients = null;        

        public List<double> time = null;
        public List<double> intensityCount = null;

        private MathNet.Numerics.Interpolation.IInterpolation interpole = null;
        
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

        public double InterpolateIntensity(double timePoint)
        {
            if (time != null && time.Count > 2)
            {
                if(interpole == null)
                    interpole = MathNet.Numerics.Interpolation.Interpolate.LinearBetweenPoints(time, intensityCount);
            
                return interpole.Interpolate(timePoint);
            }
            return 0;
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
