/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Proteomics.Utilities;

namespace Trinity
{
    public class PeptideSpectrumMatches : GraphML_List<PeptideSpectrumMatch>
    {
        FDRizer<PeptideSpectrumMatch> uptimizer;
        public PeptideSpectrumMatches() { }
        public PeptideSpectrumMatches(IEnumerable<PeptideSpectrumMatch> list) : base(list) { }

        public void OptimizePSMScoreRatios(DBOptions options, double desired_fdr, Result results)
        {
            long bestNbTargets = 0;
            double bestProtein = 0.1;
            double bestPeptide = 0.2;
            double bestFragments = 0.4;
            double bestIntensities = 0.1;
            double bestPrecursor = 0.2;

            double incr = 0.05;
            for (options.dProtein = 0.1; options.dProtein < 0.3; options.dProtein += incr)
                for (options.dPeptideScore = 0.1; options.dPeptideScore < 0.4; options.dPeptideScore += incr)
                    for (options.dPrecursor = 0.1; options.dPrecursor < 0.5; options.dPrecursor += incr)
                        for (options.dMatchingProductFraction = 0.1; options.dMatchingProductFraction < 1 - options.dPrecursor - options.dPeptideScore - options.dProtein; options.dMatchingProductFraction += incr)
                        {
                            double cumul = options.dPrecursor + options.dPeptideScore + options.dProtein + options.dMatchingProductFraction;
                            for (options.dIntensityFraction = 0.1; options.dIntensityFraction < 1 - cumul; options.dIntensityFraction += incr)
                            {
                                if (cumul + options.dIntensityFraction == 1)
                                {
                                    long nbTargets = 0;
                                    foreach (Precursor precursor in results.matchedPrecursors)
                                        if (precursor.Target)
                                            nbTargets++;
                                    if (nbTargets > bestNbTargets)
                                    {
                                        bestNbTargets = nbTargets;
                                        bestProtein = options.dProtein;
                                        bestPeptide = options.dPeptideScore;
                                        bestPrecursor = options.dPrecursor;
                                        bestIntensities = options.dIntensityFraction;
                                        bestFragments = options.dMatchingProductFraction;
                                    }
                                }
                            }

                        }

            options.dProtein = bestProtein;
            options.dPeptideScore = bestPeptide;
            options.dPrecursor = bestPrecursor;
            options.dIntensityFraction = bestIntensities;
            options.dMatchingProductFraction = bestFragments;

            options.ConSole.WriteLine("New score ratios   ----------------------------------------------------------------------- ");
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dPrecursor:                     " + options.dPrecursor);
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dMatchingProductFraction:       " + options.dMatchingProductFraction);
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dIntensityFraction:             " + options.dIntensityFraction);
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dPeptideScore:                  " + options.dPeptideScore);
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dProtein:                       " + options.dProtein);
            options.ConSole.WriteLine("------------------------------------------------------------------------------------------ ");
        }

        public void OptimizePSMScoreRatios_PREVIOUSVersion(DBOptions options, double desired_fdr)
        {
            //TODO find a Max flow approach to modelize the score optimization routine to maximize Targets versus Decoys
            List<PeptideSpectrumMatch> sortedPrecursorPrecision = new List<PeptideSpectrumMatch>(this);
            sortedPrecursorPrecision.Sort(PeptideSpectrumMatches.ComparePrecursorScore);
            double ratioPrecursorPrecision = FDRizer<PeptideSpectrumMatch>.ComputeAtFDR(sortedPrecursorPrecision, desired_fdr).Count / (double)this.Count;

            List<PeptideSpectrumMatch> sortedFragments = new List<PeptideSpectrumMatch>(this);
            sortedFragments.Sort(PeptideSpectrumMatches.CompareMatchingProductsFraction);
            double ratioFragments = FDRizer<PeptideSpectrumMatch>.ComputeAtFDR(sortedFragments, desired_fdr).Count / (double)this.Count;

            List<PeptideSpectrumMatch> sortedIntensities = new List<PeptideSpectrumMatch>(this);
            sortedIntensities.Sort(PeptideSpectrumMatches.CompareMatchingIntensityFraction);
            double ratioIntensities = FDRizer<PeptideSpectrumMatch>.ComputeAtFDR(sortedIntensities, desired_fdr).Count / (double)this.Count;

            List<PeptideSpectrumMatch> sortedPeptides = new List<PeptideSpectrumMatch>(this);
            sortedPeptides.Sort(PeptideSpectrumMatches.ComparePeptideScore);
            double ratioPeptideScore = FDRizer<PeptideSpectrumMatch>.ComputeAtFDR(sortedPeptides, desired_fdr).Count / (double)this.Count;

            double totalRatios = ratioPrecursorPrecision + ratioFragments + ratioIntensities + ratioPeptideScore;
            options.dPrecursor += ratioPrecursorPrecision / totalRatios;
            options.dMatchingProductFraction += ratioFragments / totalRatios;
            options.dIntensityFraction += ratioIntensities / totalRatios;
            options.dPeptideScore += ratioPeptideScore / totalRatios;

            double somme = options.dPrecursor + options.dMatchingProductFraction + options.dIntensityFraction + options.dPeptideScore;
            options.dPrecursor /= somme;
            options.dMatchingProductFraction /= somme;
            options.dIntensityFraction /= somme;
            options.dPeptideScore /= somme;

            //options.dPrecursor                 = ratioPrecursorPrecision   / totalRatios;
            //options.dMatchingProductFraction   = ratioFragments            / totalRatios;
            //options.dIntensityFraction         = ratioIntensities          / totalRatios;
            //options.dPeptideScore              = ratioPeptideScore         / totalRatios;
            options.ConSole.WriteLine("New score ratios  [" + totalRatios + " total ratios] ------------------------------------- ");
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dPrecursor:                     " + options.dPrecursor);
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dMatchingProductFraction:       " + options.dMatchingProductFraction);
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dIntensityFraction:             " + options.dIntensityFraction);
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dPeptideScore:                  " + options.dPeptideScore);
            options.ConSole.WriteLine("------------------------------------------------------------------------------------------ ");
        }

        public double ComputePrecursorArea(bool smooth)
        {
            double lastTimeStamp = 0;
            double lastIntensityPerUnitOfTime = 0;
            List<double> timeGap = new List<double>();
            List<double> precursorIntensities = new List<double>();
            this.Sort(PeptideSpectrumMatches.AscendingRetentionTime);

            foreach (PeptideSpectrumMatch psm in this)
            {
                if (psm.Query.spectrum.Ms1InjectionTime > 0)
                {
                    if (psm.Query.spectrum.PrecursorIntensity > 0 && lastTimeStamp > 0)
                    {
                        timeGap.Add((psm.Query.spectrum.RetentionTimeInMin - lastTimeStamp) * 60.0 * 1000.0);
                        precursorIntensities.Add(psm.Query.spectrum.PrecursorIntensityPerMilliSecond);
                    }
                    lastIntensityPerUnitOfTime = psm.Query.spectrum.PrecursorIntensityPerMilliSecond;
                    lastTimeStamp = psm.Query.spectrum.RetentionTimeInMin;
                }
            }
            //Smooth the curve

            if (smooth)
            {
                List<double> timeGapSmooth = new List<double>();
                List<double> precursorIntensitiesSmoothed = new List<double>();
                for (int i = 0; i < timeGap.Count; i++)
                {
                    double cumulIntensity = 0;
                    int nbItems = 0;
                    for (int k = (i - 2 < 0 ? 0 : i - 2); k < timeGap.Count - 2; k++)
                    {
                        nbItems++;
                        cumulIntensity += precursorIntensities[k];
                    }
                    if (nbItems > 0)
                    {
                        timeGapSmooth.Add(timeGap[0]);
                        precursorIntensitiesSmoothed.Add(cumulIntensity / (double)nbItems);
                    }
                }

                double fragSpectrumArea = 0;
                for (int i = 0; i < timeGapSmooth.Count; i++)
                    fragSpectrumArea += timeGapSmooth[i] * precursorIntensitiesSmoothed[i];
                return fragSpectrumArea;
            }
            else
            {
                double fragSpectrumArea = 0;
                for (int i = 0; i < timeGap.Count; i++)
                    fragSpectrumArea += timeGap[i] * precursorIntensities[i];
                return fragSpectrumArea;
            }
        }

        public double ComputePrecursorAreaBKP(bool smooth)
        {
            double lastTimeStamp = 0;
            double lastIntensityPerUnitOfTime = 0;
            List<double> timeGap = new List<double>();
            List<double> precursorIntensities = new List<double>();
            foreach (PeptideSpectrumMatch psm in this)
            {
                if (psm.Query.spectrum.Ms1InjectionTime > 0)
                {
                    if (psm.Query.spectrum.PrecursorIntensity > 0 && lastTimeStamp > 0)
                    {
                        timeGap.Add((psm.Query.spectrum.RetentionTimeInMin - lastTimeStamp) * 60.0 * 1000.0);
                        if (lastIntensityPerUnitOfTime > 0)
                            precursorIntensities.Add((psm.Query.spectrum.PrecursorIntensity / psm.Query.spectrum.Ms1InjectionTime + lastIntensityPerUnitOfTime) * 0.5);
                        else
                            precursorIntensities.Add(psm.Query.spectrum.PrecursorIntensity / psm.Query.spectrum.Ms1InjectionTime);
                    }
                    lastIntensityPerUnitOfTime = psm.Query.spectrum.PrecursorIntensity / psm.Query.spectrum.Ms1InjectionTime;
                    lastTimeStamp = psm.Query.spectrum.RetentionTimeInMin;
                }
            }
            //Smooth the curve

            if (smooth)
            {
                List<double> timeGapSmooth = new List<double>();
                List<double> precursorIntensitiesSmoothed = new List<double>();
                for (int i = 0; i < timeGap.Count; i++)
                {
                    double cumulIntensity = 0;
                    int nbItems = 0;
                    for (int k = (i - 2 < 0 ? 0 : i - 2); k < timeGap.Count - 2; k++)
                    {
                        nbItems++;
                        cumulIntensity += precursorIntensities[k];
                    }
                    if (nbItems > 0)
                    {
                        timeGapSmooth.Add(timeGap[0]);
                        precursorIntensitiesSmoothed.Add(cumulIntensity / (double)nbItems);
                    }
                }

                double fragSpectrumArea = 0;
                for (int i = 0; i < timeGapSmooth.Count; i++)
                    fragSpectrumArea += timeGapSmooth[i] * precursorIntensitiesSmoothed[i];
                return fragSpectrumArea;
            }
            else
            {
                double fragSpectrumArea = 0;
                for (int i = 0; i < timeGap.Count; i++)
                    fragSpectrumArea += timeGap[i] * precursorIntensities[i];
                return fragSpectrumArea;
            }
        }

        public Dictionary<PeptideSpectrumMatch, double> ComputeMsMsNormalizationFactors()
        {
            Dictionary<PeptideSpectrumMatch, double> fragRatio = new Dictionary<PeptideSpectrumMatch, double>();
            double lastIntensity = 0;
            foreach (PeptideSpectrumMatch psm in this)
            {
                if (psm.Query.spectrum.PrecursorIntensity > 0)
                {
                    double intensityFactor = 0;
                    if (psm.Query.spectrum.InjectionTime >= 119.999997317791)//TODO Find instrument/method specific default or max injection time
                        lastIntensity = psm.Query.spectrum.PrecursorIntensity;
                    else
                    {
                        double predictedIntensity = (lastIntensity);// + psm.Query.spectrum.PrecursorIntensity) * 0.5;
                        intensityFactor = (psm.Query.spectrum.PrecursorIntensity - predictedIntensity) / predictedIntensity;// psm.Query.spectrum.PrecursorIntensity;
                    }
                    fragRatio.Add(psm, intensityFactor);
                }
            }
            return fragRatio;
        }        

        public List<ProductMatch> GetCombinedSpectrum(DBOptions dbOptions, Peptide peptide, int psmCharge, Dictionary<double, int> DicOfCommonPM = null)
        {
            //Dictionary<PeptideSpectrumMatch, double> DicOfPsmFactor = this.ComputeMsMsNormalizationFactors();
            //Dictionary<ProductMatch, double> DicOfProductMsMsFactor = new Dictionary<ProductMatch, double>();
            Dictionary<string, List<ProductMatch>> DicOfProducts = new Dictionary<string, List<ProductMatch>>();

            foreach (PeptideSpectrumMatch psm in this)
            {
                foreach (ProductMatch match in psm.AllProductMatches)
                {
                    if (DicOfCommonPM == null || DicOfCommonPM.ContainsKey(match.theoMz))
                    {
                        string key = match.fragment + "|" + match.fragmentPos + "|" + match.charge;
                        if (!DicOfProducts.ContainsKey(key))
                            DicOfProducts.Add(key, new List<ProductMatch>());
                        DicOfProducts[key].Add(match);
                        //DicOfProductMsMsFactor.Add(match, DicOfPsmFactor[psm]);
                    }
                }
            }

            double avgInt = 0;
            List<ProductMatch> products = new List<ProductMatch>();
            //if (DicOfProductMsMsFactor.Count > 0)
            {
                foreach (List<ProductMatch> matchList in DicOfProducts.Values)
                {
                    ProductMatch newPM = new ProductMatch(matchList[0]);
                    newPM.obsIntensity = 0;
                    if (matchList.Count > 0)
                    {
                        foreach (ProductMatch pm in matchList)
                        {
                            newPM.obsIntensity += pm.obsIntensity;// +pm.obsIntensity * DicOfProductMsMsFactor[pm];
                            //newPM.obsIntensity += pm.obsIntensity + pm.obsIntensity * DicOfProductMsMsFactor[pm];
                            newPM.normalizedIntensity += pm.normalizedIntensity;
                        }

                        newPM.obsIntensity /= (double)matchList.Count;
                        newPM.normalizedIntensity /= (double)matchList.Count;
                    }
                    newPM.weight = matchList.Count * newPM.obsIntensity;
                    avgInt += newPM.normalizedIntensity;
                    products.Add(newPM);
                }
                avgInt /= (double)products.Count;

                //Add missed important fragments
                if (DicOfCommonPM != null)
                {
                    foreach (double mz in DicOfCommonPM.Keys)
                    {
                        bool found = false;
                        foreach (ProductMatch match in products)
                            if (match.theoMz == mz)
                                found = true;
                        if (!found)
                        {
                            ProductMatch newMatch = new ProductMatch();
                            newMatch.theoMz = mz;
                            newMatch.weight = 0;
                            newMatch.obsIntensity = 0;
                            newMatch.normalizedIntensity = 0;
                            foreach (PeptideSpectrumMatch psm in this)
                            {
                                foreach (MsMsPeak peak in psm.Query.spectrum.Peaks)
                                {
                                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(peak.MZ, mz, dbOptions.productMassTolerance.Units)) <= dbOptions.productMassTolerance.Value)
                                    {
                                        newMatch.weight += 1;
                                        newMatch.obsIntensity += peak.Intensity;// + peak.Intensity * DicOfPsmFactor[psm];
                                        newMatch.normalizedIntensity += peak.Intensity / (psm.Query.spectrum.PrecursorIntensityPerMilliSecond * psm.Query.spectrum.InjectionTime);
                                    }
                                }
                            }
                            if (newMatch.normalizedIntensity < avgInt * 0.05)
                            {
                                newMatch.normalizedIntensity = 0;
                                newMatch.obsIntensity = 0;
                            }
                            else
                            {
                                newMatch.obsIntensity /= (double)newMatch.weight;
                                newMatch.normalizedIntensity /= (double) newMatch.weight;
                            }
                            newMatch.weight *= newMatch.obsIntensity;
                            products.Add(newMatch);
                        }
                    }
                }
            }
            return products;
        }

        public static int CompareMatchingIntensityFraction(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingIntensityFraction.CompareTo(right.MatchingIntensityFraction);
        }
        public List<PeptideSpectrumMatch> ComputeAtFDR(double desired_fdr, bool precursorErrorOnly)
        {
            if (uptimizer == null)
            {
                List<Comparison<PeptideSpectrumMatch>> sorts = new List<Comparison<PeptideSpectrumMatch>>();
                if (precursorErrorOnly)
                    sorts.Add(CompareProductScore);
                else
                {
                    //sorts.Add(CompareMatchingIntensity);
                    sorts.Add(CompareMaxQuantScore);
                    //                sorts.Add(ComparePrecursorScore);
                    //                sorts.Add(CompareMatchingProductsFraction);
                    //                sorts.Add(CompareMatchingProducts);
                    //sorts.Add(CompareProteinScore);
                    //sorts.Add(ComparePeptideScore);
                }
                uptimizer = new FDRizer<PeptideSpectrumMatch>(this, sorts, null);
            }
            return uptimizer.Launch(desired_fdr);
        }
        /*
            return              dIntensity * MatchingIntensityFraction +
                                dProduct * ProductScore +
                                dPrecursor * PrecursorScore +
                                dMatchingProductFraction * MatchingProductsFraction +
                                dMatchingProduct * MatchingProducts +
                                dProtein * ProteinScore +
                                dPeptideScore * PeptideScore;*/
        public static int CompareTrapDistance(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.Query.TrapDistance.CompareTo(right.Query.TrapDistance);
        }
        public static int AscendingRetentionTime(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return left.Query.spectrum.RetentionTimeInMin.CompareTo(right.Query.spectrum.RetentionTimeInMin);
        }
        public static int CompareMatchingIntensity(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingIntensity.CompareTo(right.MatchingIntensity);
        }
        public static int CompareMaxQuantScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MaxQuantScore().CompareTo(right.MaxQuantScore());
        }
        public static int CompareProductScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.ProductScore.CompareTo(right.ProductScore);
        }
        public static int ComparePrecursorScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.PrecursorScore.CompareTo(right.PrecursorScore);
        }
        public static int CompareMatchingProductsFraction(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingProductsFraction.CompareTo(right.MatchingProductsFraction);
        }
        public static int CompareMatchingProducts(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingWeightedProducts.CompareTo(right.MatchingWeightedProducts);
        }
        public static int CompareProteinScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.ProteinScore.CompareTo(right.ProteinScore);
        }
        public static int ComparePeptideScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.PeptideScore.CompareTo(right.PeptideScore);
        }
    }
}
