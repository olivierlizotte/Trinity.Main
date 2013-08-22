/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Text;
using System.Collections.Generic;
using Proteomics.Utilities;
using Proteomics.Utilities.MaxQuant;

namespace Trinity
{
    public class ProductMatch : GraphML_Node
    {
        public double theoMz;
        public double obsMz;
        public int    charge;
        public double obsIntensity;
        public double mass_diff;
        public string fragment;
        public int    fragmentPos;
        public double weight;
        public ProductMatch()
        {
        }

        public ProductMatch(ProductMatch match)
        {
            this.theoMz = match.theoMz;
            this.obsMz = match.obsMz;
            this.obsIntensity = match.obsIntensity;
            this.mass_diff = match.mass_diff;
            this.charge = match.charge;
            this.fragment = match.fragment;
            this.fragmentPos = match.fragmentPos;
            this.weight = match.weight;
        }
    }
    public class PeptideSpectrumMatch : GraphML_Node, ITargetDecoy
    {
        public double PeptideScore;
        public double ProteinScore;
        
        public Query Query { get; set; }

        public Peptide Peptide { get; set; }

        public double PrecursorMzError { get; set; }

        public int MatchingProducts { get; set; }
        public double MatchingWeightedProducts { get; set; }
        public double MatchingProductsProbability { get; set; }
        public GraphML_List<ProductMatch> AllProductMatches{ get; set; }

        public int TotalTheoreticalProducts { get; set; }
        public double TotalWeightedProducts { get; set; }

        public double MatchingProductsFraction { get; set; }

        public double MatchingIntensity { get; set; }

        public double MatchingIntensityFraction { get; set; }

        public double MorpheusScore { get;  set; }
        public double PrecursorScore { get;  set; }
        public double ProductScore { get;  set; }
        public double IntensityScore { get;  set; }

        //TODO Check if those var are common with Precursor and remove zeroed out var  
        public double ProbabilityScore()
        {
            DBOptions options = Query.options;
            //return MaxQuantScore();//TODO Reactivate Optimized PeptideSpectrumMatch Score
            double score = options.dIntensity * MatchingIntensity +
                                options.dIntensityFraction * MatchingIntensityFraction +
                                options.dProduct * ProductScore +
                                options.dPrecursor * PrecursorScore +
                                options.dMatchingProductFraction * MatchingProductsFraction +
                                options.dMatchingProduct * MatchingWeightedProducts +
                                options.dProtein * ProteinScore +
                                options.dPeptideScore * PeptideScore;//*/
            return score;
        }

        public bool Decoy
        {
            get { return Peptide.Decoy; }
        }

        public bool Target
        {
            get { return !Decoy; }
        }

        private static readonly ProductTypes PRODUCT_TYPES = ProductTypes.Instance;

        public PeptideSpectrumMatch()
        {
            AllProductMatches = new GraphML_List<ProductMatch>();
        }

        public PeptideSpectrumMatch(Query query, Peptide peptide, DBOptions options)
        {
            Query = query;
            Peptide = peptide;

            UpdatePrecursor(options);

            Initialize(options, GetProductMZs(options, query.spectrum.Peaks), query.spectrum.MostIntensePeak);
        }

        public void UpdatePrecursor(DBOptions options)
        {
            //PrecursorMassErrorDa = (spectrum. - spectrum.IsotopicShift) - peptide.MonoisotopicMass;
            PrecursorMzError = Numerics.MzDifference(Numerics.MZFromMass(Peptide.MonoisotopicMass, Query.precursor.Charge), Query.precursor.Track.MZ, options.precursorMassTolerance.Units);
            if (PrecursorMzError > options.precursorMassTolerance.Value)
                PrecursorMzError = options.precursorMassTolerance.Value;

            PrecursorScore = (options.precursorMassTolerance.Value - Math.Abs(PrecursorMzError)) / options.precursorMassTolerance.Value;
        }
        
        public void Merge(PeptideSpectrumMatch psm, DBOptions options)
        {
            //Augment list A
            List<ProductMatch> newList = new List<ProductMatch>();
            foreach (ProductMatch matchA in this.AllProductMatches)
            {
                bool isNew = true;
                foreach (ProductMatch matchB in psm.AllProductMatches)
                    if (matchA.theoMz == matchB.theoMz && matchA.charge == matchB.charge)
                    {
                        if (matchA.mass_diff > matchB.mass_diff)
                        {
                            matchA.mass_diff = matchB.mass_diff;
                            matchA.obsMz = matchB.obsMz;
                        }
                        if (matchA.obsIntensity < matchB.obsIntensity)
                            matchA.obsIntensity = matchB.obsIntensity;
                        newList.Add(matchA);
                        isNew = false;
                    }
                if (isNew)
                    newList.Add(matchA);
            }
            //Add items unique to list B
            foreach (ProductMatch matchB in psm.AllProductMatches)
            {
                bool isNew = true;
                foreach (ProductMatch matchNew in newList)
                    if (matchNew.theoMz == matchB.theoMz && matchNew.charge == matchB.charge)
                        isNew = false;
                if (isNew)
                    newList.Add(matchB);
            }
            double mostIntenseFragment = this.Query.spectrum.MostIntensePeak;
            if (mostIntenseFragment < psm.Query.spectrum.MostIntensePeak)
                mostIntenseFragment = psm.Query.spectrum.MostIntensePeak;
            Initialize(options, newList, mostIntenseFragment);
        }//*/

        public IEnumerable<ProductMatch> GetProductMZs(DBOptions options, GraphML_List<MsMsPeak> peaks)//, List<double> theoretical_product_mzs = null)
        {
            //if (theoretical_product_mzs == null)
                //theoretical_product_mzs = Peptide.CalculateAllProductMz(PRODUCT_TYPES[Query.spectrum.FragmentationMethod], Query.precursor);

            // speed optimizations
            //double[] experimental_masses = Query.spectrum.Masses;
            //GraphML_List<MSPeak> peaks = Query.spectrum.Peaks;
            int num_experimental_peaks = peaks.Count;
            double product_mass_tolerance_value = options.productMassTolerance.Value;
            MassToleranceUnits product_mass_tolerance_units = options.productMassTolerance.Units;
            TotalTheoreticalProducts = 0;
            TotalWeightedProducts = 0;
            //New version that should include charged ions
            ProductMatch pMatch = new ProductMatch();
            foreach (ProductMatch matchTheo in options.fragments.ComputeFragments(Peptide, Query.precursor.Charge, options))//Peptide.EnumerateAllProductMz(PRODUCT_TYPES[Query.spectrum.FragmentationMethod], Query.precursor))
            //foreach (ProductMatch matchTheo in Peptide.EnumerateAllProductMz(PRODUCT_TYPES[Query.spectrum.FragmentationMethod], Query.precursor))
            {
                TotalTheoreticalProducts++;
                TotalWeightedProducts += matchTheo.weight;//++;
                //TODO test what is most common between charged ions, and uncharged ions
                //for (int charge = 1; charge <= Query.precursor.Charge; charge++)
                //for (int charge = Query.precursor.Charge; charge >= 1; charge--)
                //{
                double massDiff = options.productMassTolerance.Value;
                double bestMz = -1;
                double bestInt = -1;

                foreach (int index in Query.spectrum.GetIndexOfMZInRange(matchTheo.theoMz, options.productMassTolerance))
                {
                    if (peaks[index].Charge == 0 || peaks[index].Charge == matchTheo.charge)
                    {
                        double diff = matchTheo.theoMz - peaks[index].MZ;// experimental_masses[index];
                        if (Math.Abs(diff) < options.productMassTolerance.Value)
                        {
                            if (Math.Abs(diff) < Math.Abs(massDiff))//TODO Priority to intensity, or precision?
                            {
                                massDiff = diff;
                                bestMz = peaks[index].MZ;
                            }
                            if (peaks[index].Intensity > bestInt)
                                bestInt = peaks[index].Intensity;
                        }
                    }
                }
                if (bestMz >= 0)
                {
                    pMatch.weight = matchTheo.weight;
                    pMatch.theoMz = matchTheo.theoMz;// Utilities.MZFromMzSingleCharge(theoMass, charge);
                    pMatch.obsMz = bestMz;// experimental_masses[bestIndex];
                    pMatch.mass_diff = massDiff;
                    pMatch.obsIntensity = bestInt;// Intensities[bestIndex];
                    pMatch.charge = matchTheo.charge;// peaks[bestIndex].Charge;
                    pMatch.fragment = matchTheo.fragment;
                    pMatch.fragmentPos = matchTheo.fragmentPos;
                    yield return pMatch;
                    //break;
                }
            }
        }
        
        public double MaxQuantScore()
        {
			//Test which is better between without and with neutral loss ... keep best score
			double lnp = Math.Log(0.06);
			double lnq = Math.Log(1.0 - 0.06);
            int n = TotalTheoreticalProducts;
            int k = MatchingProducts;
			n = Math.Max(n, k);
			double s = -k * lnp - (n - k) * lnq -
                       NumericalRecipes.Factln(n) + NumericalRecipes.Factln(k) + NumericalRecipes.Factln(n - k);
            s = 10.0 * s / Math.Log(10);
			return s;		
        }

        public void Initialize(DBOptions options, IEnumerable<ProductMatch> productMZs, double highestIntensity)
        {
            //List<double> theoretical_product_mz = Peptide.CalculateAllProductMz(PRODUCT_TYPES[Query.spectrum.FragmentationMethod], Query.precursor);
            //TotalProducts = theoretical_product_mz.Count;
            MatchingProducts    = 0;
            MatchingIntensity   = 0.0;
            double cumulDiff    = 0;
            ProductScore        = 0;
            AllProductMatches = new GraphML_List<ProductMatch>();
            MatchingWeightedProducts = 0;
            foreach (ProductMatch match in productMZs)
            {
                MatchingProducts++;
                MatchingWeightedProducts  += match.weight;
                MatchingIntensity += match.obsIntensity;
                cumulDiff += Math.Abs(match.mass_diff) * match.weight;
                AllProductMatches.Add(new ProductMatch(match));
            }
            MatchingProductsFraction    = (double)MatchingWeightedProducts / (double) TotalWeightedProducts;
            MatchingIntensityFraction   = MatchingIntensity / (double)(highestIntensity * TotalTheoreticalProducts);
            if (MatchingIntensityFraction > 1)
                MatchingIntensityFraction = 1.0;
            ProductScore                = 1.0 - (cumulDiff / (double)(MatchingProducts * options.productMassTolerance.Value));
            IntensityScore              = MatchingIntensityFraction / (double)Query.spectrum.Peaks.Count;
            MorpheusScore               = MatchingProducts + MatchingIntensityFraction;
        }
        
        public static int AscendingSpectrumNumberComparison(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return left.Query.spectrum.ScanNumber.CompareTo(right.Query.spectrum.ScanNumber);
        }

        public static int DescendingOptimizedScoreComparison(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            int comparison = -(left.ProbabilityScore().CompareTo(right.ProbabilityScore()));
            if (comparison != 0)
            {
                return comparison;
            }
            else
            {
                return left.Target.CompareTo(right.Target);
            }
        }

        public static readonly string Header = "Filename\tScan Number\tRetention Time (min)\tPrecursor m/z\tPrecursor Intensity\tPrecursor Charge\tPrecursor Mass (Da)\tExperimental Peaks\tTotal Intensity"
            + "\tPeptide Sequence\tBase Peptide Sequence\tProtein Description\tStart Residue Number\tStop Residue Number\tMissed Cleavages"
            + "\tTheoretical Mass (Da)\tPrecursor Mass Error (Da)\tPrecursor Mass Error (ppm)"
            + "\tMatching Products\tTotal Products\tRatio of Matching Products\tMatching Intensity\tFraction of Intensity Matching\tMorpheus Score";

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            //sb.Append(Query.spectrum.Filename + '\t');
            sb.Append(Query.spectrum.ScanNumber.ToString() + '\t');
            sb.Append(Query.spectrum.RetentionTimeInMin.ToString() + '\t');
            sb.Append(Query.spectrum.PrecursorMZ.ToString() + '\t');
            sb.Append(Query.spectrum.PrecursorIntensity.ToString() + '\t');
            sb.Append(Query.spectrum.PrecursorCharge.ToString() + '\t');
            sb.Append(Query.spectrum.PrecursorMass.ToString() + '\t');
            sb.Append(Query.spectrum.Peaks.Count.ToString() + '\t');
            sb.Append(Query.spectrum.TotalIntensity.ToString() + '\t');
            sb.Append(Peptide.ExtendedSequence() + '\t');
            sb.Append(Peptide.BaseSequence + '\t');
            sb.Append(Peptide.Parent.Description + '\t');
            sb.Append(Peptide.StartResidueNumber.ToString() + '\t');
            sb.Append(Peptide.EndResidueNumber.ToString() + '\t');
            sb.Append(Peptide.MissedCleavages.ToString() + '\t');
            sb.Append(Peptide.MonoisotopicMass.ToString() + '\t');
            sb.Append(PrecursorMzError.ToString() + '\t');
            sb.Append(MatchingProducts.ToString() + '\t');
            sb.Append(TotalTheoreticalProducts.ToString() + '\t');
            sb.Append(MatchingProductsFraction.ToString() + '\t');
            sb.Append(MatchingIntensity.ToString() + '\t');
            sb.Append(MatchingIntensityFraction.ToString() + '\t');
            sb.Append(ProbabilityScore().ToString());

            return sb.ToString();
        }
    }
}