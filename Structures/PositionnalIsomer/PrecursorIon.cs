using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Trinity.Structures.PositionnalIsomer
{
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
            this.eCurve = ElutionCurve.Create(this.Queries);
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
            return DicOfSpectrumMasses;
        }
    }
}
