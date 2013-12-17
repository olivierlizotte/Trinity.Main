using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Trinity.Methods;

namespace Trinity.Structures.PositionnalIsomer
{
    public class CharacterizedPrecursor : PrecursorIon
    {
        public Peptide Peptide;
        public PeptideSpectrumMatches Psms;
        public Dictionary<int, double> PrecursorLossNormalizeFactor;        
        public Dictionary<int, List<ProductMatch>> Fragments;
        public List<ProductMatch> AllFragments;

        private Dictionary<PeptideSpectrumMatch, double> DicOfPsmFactor;

        public Dictionary<int, Dictionary<double, double>> NormalizedFragments;

        public CharacterizedPrecursor(Sample sample, DBOptions dbOptions, Peptide peptide, IEnumerable<Query> queries, double mz)
            : base(sample, queries, mz, -1)
        {
            this.Peptide = peptide;
            Psms = new PeptideSpectrumMatches();

            foreach (Query query in queries)
                if (query.sample == sample)
                {
                    Psms.Add(new PeptideSpectrumMatch(query, peptide, dbOptions));
                    Charge = query.precursor.Charge;
                }
            Psms.Sort(PeptideSpectrumMatches.AscendingRetentionTime);
            AllFragments = Psms.GetCombinedSpectrum(dbOptions);
            AllFragments.Sort(ProductMatch.DescendingWeightComparison);

            DicOfPsmFactor = new Dictionary<PeptideSpectrumMatch, double>();
            foreach (PeptideSpectrumMatch psm in Psms)
                DicOfPsmFactor.Add(psm, psm.ProbabilityScore());

            Fragments = new Dictionary<int, List<ProductMatch>>();
            PrecursorLossNormalizeFactor = new Dictionary<int, double>();
        }

        private static Dictionary<double, int> GetCommonFragmentMz(IEnumerable<CharacterizedPrecursor> precursors, int nbFragmentPerPrec)
        {
            Dictionary<double, int> DicOfFragments = new Dictionary<double, int>();
            foreach (CharacterizedPrecursor cPrec in precursors)
            {
                for (int i = 0; i < nbFragmentPerPrec && i < cPrec.AllFragments.Count; i++)
                    if (!DicOfFragments.ContainsKey(cPrec.AllFragments[i].theoMz))
                        DicOfFragments.Add(cPrec.AllFragments[i].theoMz, 1);
                    else
                        DicOfFragments[cPrec.AllFragments[i].theoMz]++;
            }
            return DicOfFragments;
        }

        private List<ProductMatch> GetCombinedMatches(Dictionary<double, int> dicOfCommonFragments, DBOptions dbOptions)
        {
            List<ProductMatch> matches = new List<ProductMatch>(dicOfCommonFragments.Count);

            foreach (double mz in dicOfCommonFragments.Keys)
            {
                bool found = false;
                foreach (ProductMatch match in AllFragments)
                    if (match.theoMz == mz)
                    {
                        matches.Add(new ProductMatch(match));
                        found = true;
                    }

                if (!found)
                {
                    double sumPsmFactor = 0;
                    ProductMatch newMatch = new ProductMatch();
                    newMatch.theoMz = mz;
                    newMatch.weight = 0;
                    newMatch.obsIntensity = 0;
                    newMatch.normalizedIntensity = 0;
                    foreach (PeptideSpectrumMatch psm in Psms)
                    {
                        foreach (MsMsPeak peak in psm.Query.spectrum.Peaks)
                        {
                            if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(peak.MZ, mz, dbOptions.productMassTolerance.Units)) <= dbOptions.productMassTolerance.Value)
                            {
                                newMatch.obsIntensity += peak.Intensity * DicOfPsmFactor[psm];
                                newMatch.normalizedIntensity += (peak.Intensity / (psm.Query.spectrum.PrecursorIntensityPerMilliSecond * psm.Query.spectrum.InjectionTime)) * DicOfPsmFactor[psm];
                                sumPsmFactor += DicOfPsmFactor[psm];
                                newMatch.weight++;
                            }
                        }
                    }
                    if (newMatch.weight > 0)
                    {
                        newMatch.obsIntensity /= sumPsmFactor;
                        newMatch.normalizedIntensity /= sumPsmFactor;
                    }
                    newMatch.weight *= newMatch.normalizedIntensity;
                    matches.Add(newMatch);
                }
            }

            double averageNormedIntensity = 0.0;
            foreach (ProductMatch match in matches)
                averageNormedIntensity += match.normalizedIntensity;

            if (matches.Count > 0)
                averageNormedIntensity /= (double)matches.Count;

            //Keep only most intense fragments (5% of average normalized intensity)
            foreach (ProductMatch pm in matches)
                if (pm.normalizedIntensity < averageNormedIntensity * 0.1)//0.05
                {
                    pm.normalizedIntensity = 0;
                    pm.obsIntensity = 0;
                }

            return matches;
        }

        private double GetNormalizePrecursorFactor(IEnumerable<CharacterizedPrecursor> allCorrespondingPrec, out double average, out bool keep)
        {
            keep = false;
            double PrecursorLossNormalizeFactor = 1.0;
            average = eCurve.Area;

            //average = area * Norm => Norm = average/area
            if (eCurve.Area > 0)
            {
                //Normalize matches based on precursor intensity differences
                double cumulArea = 0.0;
                int nbNonZero = 0;
                foreach (CharacterizedPrecursor precursor in allCorrespondingPrec)
                {
                    if (precursor.eCurve.Area > 0)
                    {
                        nbNonZero++;
                        cumulArea += precursor.eCurve.Area;
                    }
                }
                if (nbNonZero > 0)
                {
                    keep = true;
                    average = cumulArea / (double)nbNonZero;
                    //PrecursorLossNormalizeFactor = Math.Log(average, 2) / Math.Log(this.eCurve.Area, 2);
                    PrecursorLossNormalizeFactor = Math.Log(average, 10) / Math.Log(this.eCurve.Area, 10);
                    //PrecursorLossNormalizeFactor = average / this.eCurve.Area;
                    //PrecursorLossNormalizeFactor = 1.0;
                }
            }

            if (PrecursorLossNormalizeFactor > 4) PrecursorLossNormalizeFactor = 4;
            if (PrecursorLossNormalizeFactor < 0.25) PrecursorLossNormalizeFactor = 0.25;

            return PrecursorLossNormalizeFactor;
        }

        private double GetNbTimesFit(DBOptions dbOptions, int nbProductsToKeep, long precision, out double area)
        {
            long nbRatios = 0;
            double nbTimes = 0;
            area = 0.0;
            List<CharacterizedPrecursor> Isomers = new List<CharacterizedPrecursor>();
            Isomers.Add(this);
            
            ElutionCurve curve = new ElutionCurve();
            foreach (Query query in this.Queries)
            {
                //double overFlow = 0;
                double underFlow = 0;
                double percentError = 0;
                Dictionary<CharacterizedPrecursor, PositionnalIsomerSolver.SolvedResult> finalRatios = PositionnalIsomerSolver.SolveFromSpectrum(Isomers, nbProductsToKeep, precision, query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                        query.spectrum.PrecursorIntensityPerMilliSecond * query.spectrum.InjectionTime, out underFlow, out percentError, dbOptions.ConSole);

                if (percentError < 0.5)
                {
                    curve.AddPoint(query.spectrum.RetentionTimeInMin * 1000.0 * 60.0, finalRatios[this].Ratio * query.spectrum.PrecursorIntensityPerMilliSecond);
                    nbTimes += finalRatios[this].NbFitTimes;
                    nbRatios++;
                }
            }
            if (nbRatios > 2 && curve.intensityCount.Count > this.Queries.Count * 0.5)
            {
                curve.Compute();
                if (curve.Area > 0)
                    area = curve.Area;
            }
            return nbTimes / (double)nbRatios;
        }

        public static void Update(IEnumerable<CharacterizedPrecursor> isomers, int minNbProducts, int maxNbProducts, DBOptions dbOptions, long precision)
        {
            foreach (CharacterizedPrecursor prec in isomers)
                prec.NormalizedFragments = new Dictionary<int,Dictionary<double,double>>();

            for (int nbProduct = minNbProducts; nbProduct <= maxNbProducts; nbProduct++)
            {
                Dictionary<double, int> dicOfFrags = GetCommonFragmentMz(isomers, nbProduct);
                foreach (CharacterizedPrecursor prec in isomers)
                    prec.Fragments.Add(nbProduct, prec.GetCombinedMatches(dicOfFrags, dbOptions));

                foreach (CharacterizedPrecursor prec in isomers)
                {
                    if (!prec.NormalizeFragments(isomers, nbProduct, dbOptions, true, false, precision))//If normalization fails, ignore this product
                        prec.Fragments.Remove(nbProduct);
                }
                
                foreach (CharacterizedPrecursor prec in isomers)
                {
                    if (prec.Fragments.ContainsKey(nbProduct))
                    {
                        Dictionary<double, double> dic = new Dictionary<double, double>();
                        foreach (double key in dicOfFrags.Keys)
                        {
                            dic.Add(key, 0.0);
                            foreach (ProductMatch match in prec.Fragments[nbProduct])
                                if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(match.theoMz, key, dbOptions.productMassTolerance.Units)) <= dbOptions.productMassTolerance.Value)
                                    dic[key] += match.normalizedIntensity;
                        }
                        prec.NormalizedFragments.Add(nbProduct, dic);
                    }
                }
            }
        }

        public static Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> GetSpikedPrecursors(Samples spikedSamples, Result spikedResult, DBOptions dbOptions, int nbMinFragments, int nbMaxFragments, long precision)
        {
            Dictionary<double, Dictionary<Query, int>> mzKeys = new Dictionary<double, Dictionary<Query, int>>();
            foreach (Query query in spikedResult.queries)
            {
                foreach (PeptideSpectrumMatch psm in query.psms)
                {
                    double mz = Proteomics.Utilities.Numerics.MZFromMass(psm.Peptide.MonoisotopicMass, query.spectrum.PrecursorCharge);
                    if (!mzKeys.ContainsKey(mz))
                    {
                        bool found = false;
                        foreach (double key in mzKeys.Keys)
                            if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(mz, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                            {
                                mz = key;
                                found = true;
                            }
                        if (!found)
                            mzKeys.Add(mz, new Dictionary<Query, int>());
                    }
                    if (!mzKeys[mz].ContainsKey(query))
                        mzKeys[mz].Add(query, 1);
                }
            }

            Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> spikes = new Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>>();
            foreach (Sample spikedSample in spikedSamples)
            {
                Dictionary<double, PrecursorIon> DicOfSpectrumMasses = PrecursorIon.GetPrecursors(spikedResult, spikedSample, dbOptions, mzKeys.Keys);
                foreach (double mzKey in DicOfSpectrumMasses.Keys)
                {
                    if (mzKeys.ContainsKey(mzKey))
                    {
                        //Pick the best PSM for each sample/precursor pair
                        Dictionary<Peptide, double> DicOfProbabilityScores = new Dictionary<Peptide, double>();

                        foreach (Query query in mzKeys[mzKey].Keys)
                            if (query.sample == spikedSample)
                            {
                                foreach (PeptideSpectrumMatch psm in query.psms)
                                    if (!DicOfProbabilityScores.ContainsKey(psm.Peptide))
                                        DicOfProbabilityScores.Add(psm.Peptide, psm.ProbabilityScore());
                                    else
                                        DicOfProbabilityScores[psm.Peptide] += psm.ProbabilityScore();
                            }

                        Peptide bestPeptide = null;
                        double bestScore = double.MinValue;
                        foreach (Peptide keyPep in DicOfProbabilityScores.Keys)
                            if (DicOfProbabilityScores[keyPep] > bestScore)
                            {
                                bestScore = DicOfProbabilityScores[keyPep];
                                bestPeptide = keyPep;
                            }
                        if (bestPeptide != null)
                        {
                            CharacterizedPrecursor cPrec = new CharacterizedPrecursor(spikedSample, dbOptions, bestPeptide, mzKeys[mzKey].Keys, mzKey);
                            //Don't keep precursors if they are not well characterized (unfragmented or missasigned)
                            if (cPrec.AllFragments.Count >= cPrec.Peptide.Length - 2)
                            {
                                if (!spikes.ContainsKey(mzKey))
                                    spikes.Add(mzKey, new Dictionary<Sample, CharacterizedPrecursor>());
                                if (!spikes[mzKey].ContainsKey(spikedSample))
                                    spikes[mzKey].Add(spikedSample, cPrec);
                                else
                                    Console.WriteLine("Twice??");
                            }
                        }
                    }
                }//End of foreach mzKey
            }//End of foreach spiked sample

            //Normalize intensities based on average area of each precursor
            List<double> tmpKeys = new List<double>(spikes.Keys);
            foreach (double mzKey in tmpKeys)
            {
                if(spikes[mzKey].Count > 1)
                    CharacterizedPrecursor.Update(spikes[mzKey].Values, nbMinFragments, nbMaxFragments, dbOptions, precision);
                else
                    spikes.Remove(mzKey);
            }//*/
            return spikes;
        }

        private bool NormalizeFragments(IEnumerable<CharacterizedPrecursor> allCorrespondingPrec, int nbProductsToKeep, DBOptions dbOptions, bool normalizePrecursor, bool normalizeFragments, long precision)
        {
            bool keepNbProds = false;
            double FragmentLossNormalizeFactor = 1.0;
            if (eCurve.Area > 0)
            {
                if (normalizeFragments)
                {
                    double area = this.eCurve.Area;
                    int nbIter = 3;
                    while (nbIter > 0)
                    {
                        nbIter--;

                        double averageNbTimes = GetNbTimesFit(dbOptions, nbProductsToKeep, precision, out area);
                        FragmentLossNormalizeFactor = precision / (double)averageNbTimes;
                        foreach (ProductMatch pm in this.Fragments[nbProductsToKeep])
                        {
                            pm.normalizedIntensity /= FragmentLossNormalizeFactor;
                            pm.obsIntensity /= FragmentLossNormalizeFactor;
                        }
                    }
                }
                else
                    FragmentLossNormalizeFactor = 1.0;

                if (normalizePrecursor)
                {
                    double average = 0;
                    PrecursorLossNormalizeFactor[nbProductsToKeep] = GetNormalizePrecursorFactor(allCorrespondingPrec, out average, out keepNbProds);
                }
                else
                    PrecursorLossNormalizeFactor[nbProductsToKeep] = 1.0;
            }
            return keepNbProds;
        }
    }
}
