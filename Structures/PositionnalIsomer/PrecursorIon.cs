using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Proteomics.Utilities.Methods;

namespace Trinity.Structures.PositionnalIsomer
{
    public class Ion
    {
        public double MZ;
        public double RetentionTime;
        public Ion(double mz, double retention)
        {
            MZ = mz;
            RetentionTime = retention;
        }
    }

    public class PrecursorIon
    {
        public ElutionCurve eCurve;
        public double MZ;
        public int Charge;
        public Queries Queries;
        public Sample Sample;
        public PrecursorIon(Sample sample, IEnumerable<Query> queries, double mz, int charge)
        {
            this.MZ = mz;
            this.Charge = charge;
            this.Queries = new Trinity.Queries();
            foreach (Query query in queries)
                if (query.sample == sample)
                    this.Queries.Add(query);

            this.Queries.Sort(Query.AscendingRetentionTimeComparison);
            Dictionary<double, double> dicOfTimeInMsVsIntensityPerMs = new Dictionary<double, double>();
            foreach (Query query in this.Queries)
                dicOfTimeInMsVsIntensityPerMs.Add(query.spectrum.RetentionTimeInMin * 60.0 * 1000.0, query.spectrum.PrecursorIntensity);//query.spectrum.PrecursorIntensityPerMilliSecond);
            this.eCurve = ElutionCurve.Create(dicOfTimeInMsVsIntensityPerMs);
            this.Sample = sample;
        }
        
        public static Dictionary<double, PrecursorIon> GetPrecursors(Result result, Sample sample, DBOptions dbOptions, IEnumerable<double> keys)
        {
            Dictionary<double, PrecursorIon> DicOfSpectrumMasses = new Dictionary<double, PrecursorIon>();            
            foreach (Query query in result.queries)
            {
                if (query.sample == sample)
                {
                    double foundKey = query.spectrum.PrecursorMZ;

                    bool foundInKeys = false;
                    foreach (double key in keys)
                    {
                        double distance = Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(query.spectrum.PrecursorMZ, key, dbOptions.precursorMassTolerance.Units));
                        if (distance <= dbOptions.precursorMassTolerance.Value)
                        {
                            if (!foundInKeys || distance < Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(query.spectrum.PrecursorMZ, foundKey, dbOptions.precursorMassTolerance.Units)))
                                foundKey = key;
                            foundInKeys = true;
                        }
                    }
                    if (!foundInKeys)
                        foreach (double key in DicOfSpectrumMasses.Keys)
                            if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(query.spectrum.PrecursorMZ, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                foundKey = key;

                    if (!DicOfSpectrumMasses.ContainsKey(foundKey))
                    {
                        PrecursorIon precIon = new PrecursorIon(sample, new Queries(dbOptions), foundKey, -1);
                        precIon.Queries.Add(query);
                        DicOfSpectrumMasses.Add(foundKey, precIon);
                    }
                    else
                        DicOfSpectrumMasses[foundKey].Queries.Add(query);
                }
            }

            //Split similar precursor mass not eluting at the same time

            //aussi://retirer le processus de clustering de propheus
            return DicOfSpectrumMasses;
        }

        public IEnumerable<PrecursorIon> SplitBasedOnTime(DBOptions dbOptions)
        {
            if(Queries.Count > 0)
            {
                this.Queries.Sort(Query.AscendingRetentionTimeComparison);
                List<double> timePoints = new List<double>();
                for (int i = 1; i < Queries.Count; i++)
                    timePoints.Add(Queries[i].spectrum.RetentionTimeInMin - Queries[i-1].spectrum.RetentionTimeInMin);

                double variance = MathNet.Numerics.Statistics.Statistics.UpperQuartile(timePoints);
                Queries newQ = new Queries(dbOptions);
                newQ.Add(Queries[0]);
                for(int i = 1; i < Queries.Count; i++)
                {
                    if(timePoints[i-1] > 10 * variance)
                    {
                        yield return new PrecursorIon(Sample, newQ, MZ, Charge);
                        newQ.Clear();
                    }
                    newQ.Add(Queries[i]);
                }
                yield return new PrecursorIon(Sample, newQ, MZ, Charge);
            }
        }
    }
}
