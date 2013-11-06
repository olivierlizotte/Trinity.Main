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
    /// <summary>
    /// Stores alignement associated functions.
    /// These method use PolynomialRegression to fit a quadratic curve for alignement of Mz values.
    /// Can be used to align Precursors (MS) or Fragments (MSMS) to fit theoretical values
    /// </summary>
    public static class Align
    {
        /// <summary>
        /// Aligns fragments observed Mz values
        /// </summary>
        /// <param name="result"></param>
        /// <param name="allPSMs"></param>
        public static void AlignProductsByDiff(Result result, PeptideSpectrumMatches allPSMs)
        {
            List<double> observedDiff = new List<double>();
            List<double> observedMz = new List<double>();
            //Precursors and Tracks
            foreach (Precursor precursor in result.matchedPrecursors)//.ComputeAtFDR(result.dbOptions.maximumFalseDiscoveryRate, false))
            {
                PeptideSpectrumMatch psm = precursor.OptimizedBestPsm();
                if (psm.Target)
                    foreach(ProductMatch fragment in psm.AllProductMatches)
                    {
                        observedMz.Add(fragment.obsMz);
                        observedDiff.Add(fragment.mass_diff);
                    }
            }
            PolynominalRegression pr = new PolynominalRegression(observedMz, observedDiff, 2);
            foreach (PeptideSpectrumMatch psm in allPSMs)
            {
                foreach (ProductMatch match in psm.AllProductMatches)
                {
                    match.obsMz += pr.Calculate(match.obsMz);
                    match.mass_diff = match.theoMz - match.obsMz;
                }
                psm.Initialize(result.dbOptions, psm.AllProductMatches);

                foreach (MsMsPeak peak in psm.Query.spectrum.Peaks)
                    peak.MZ += pr.Calculate(peak.MZ);
            }
        }
        
        /// <summary>
        /// Aligns precursors observed Mz values
        /// </summary>
        /// <param name="result"></param>
        /// <param name="allPSMs"></param>
        public static void AlignPrecursorsByDiff(Result result, PeptideSpectrumMatches allPSMs)
        {
            List<double> observedDiff = new List<double>();
            List<double> observedMz = new List<double>();
            //Precursors and Tracks
            foreach (Precursor precursor in result.matchedPrecursors)//.ComputeAtFDR(result.dbOptions.maximumFalseDiscoveryRate, false))
            {
                PeptideSpectrumMatch psm = precursor.OptimizedBestPsm();
                if (psm.Target)
                {
                    observedMz.Add(precursor.Track.MZ);
                    observedDiff.Add(Numerics.MZFromMass(psm.Peptide.MonoisotopicMass, precursor.Charge) - precursor.Track.MZ);
                }
            }
            PolynominalRegression pr = new PolynominalRegression(observedMz, observedDiff, 2);
            foreach (Query query in result.queries)
            {
                query.precursor.Track.MZ += pr.Calculate(query.precursor.Track.MZ);
                query.precursor.Mass = Numerics.MassFromMZ(query.precursor.Track.MZ, query.precursor.Charge);
                foreach (Precursor precursor in query.precursor.Isotopes)
                {
                    precursor.Track.MZ += pr.Calculate(precursor.Track.MZ);
                    precursor.Mass = Numerics.MassFromMZ(precursor.Track.MZ, precursor.Charge);
                }
                foreach (Precursor precursor in query.precursor.OtherCharges)
                {
                    precursor.Track.MZ += pr.Calculate(precursor.Track.MZ);
                    precursor.Mass = Numerics.MassFromMZ(precursor.Track.MZ, precursor.Charge);
                }

                query.spectrum.PrecursorMZ += pr.Calculate(query.spectrum.PrecursorMZ);
                query.spectrum.PrecursorMass = Numerics.MassFromMZ(query.spectrum.PrecursorMZ, query.spectrum.PrecursorCharge);
            }
            foreach (PeptideSpectrumMatch psm in allPSMs)
                psm.UpdatePrecursor(result.dbOptions);
        }

        /// <summary>
        /// Crop observed fragments outside of the variance or the standard deviation (must be below the biggest of either value)
        /// Will alter DBOptions.productMassTolerance value to reflect the change
        /// </summary>
        /// <param name="result"></param>
        /// <param name="allPSMs"></param>
        /// <returns>Returns the newly computed Fragment/Product tolerance</returns>
        public static double CropProducts(Result result, PeptideSpectrumMatches allPSMs)
        {
            List<double> errorProduct = new List<double>(result.precursors.Count);
            foreach (Precursor precursor in result.matchedPrecursors)
            {
                PeptideSpectrumMatch psm = precursor.OptimizedBestPsm();
                if (psm.Target)
                    foreach (ProductMatch match in psm.AllProductMatches)
                        errorProduct.Add(match.mass_diff);
            }
            double variance = Numerics.Variance(errorProduct);
            double stdev = Numerics.StandardDeviation(errorProduct);
            result.dbOptions.ConSole.WriteLine("Computed Product Variance = " + variance + "          STDev = " + stdev);
            if (variance < stdev)
                variance = stdev;
            //variance = result.dbOptions.productMassTolerance.Value * ((2 * variance) / result.dbOptions.productMassTolerance.Value);            

            int nbRemovedProduct = 0;
            foreach (PeptideSpectrumMatch psm in allPSMs)
            {
                for (int i = 0; i < psm.AllProductMatches.Count; )
                {
                    if (Math.Abs(psm.AllProductMatches[i].mass_diff) > variance)
                    {
                        psm.AllProductMatches.RemoveAt(i);                        
                        nbRemovedProduct++;
                    }
                    else
                        i++;
                }
                psm.MatchingProducts = psm.AllProductMatches.Count;
            }

            int nbRemovedPSM = 0;
            foreach (Precursor precursor in result.precursors)
            {
                for (int i = 0; i < precursor.psms.Count; )
                {
                    if (precursor.psms[i].MatchingProducts < 2)
                        precursor.psms.RemoveAt(i);
                    else
                    {
                        precursor.psms[i].Initialize(result.dbOptions, precursor.psms[i].AllProductMatches);
                        i++;
                    }
                }
            }            
            result.dbOptions.ConSole.WriteLine("Removed " + nbRemovedProduct + " [" + nbRemovedPSM + " removed PSMs] Fragment matches outside the variance [" + variance + "]");
            return variance;
        }

        /// <summary>
        /// Crop observed precursors values outside the variance or Standard Deviation (whichever is bigger)
        /// </summary>
        /// <param name="result"></param>
        /// <param name="allPSMs"></param>
        /// <returns>Returns the newly computed Precursor tolerance</returns>
        public static double CropPrecursors(Result result, PeptideSpectrumMatches allPSMs)
        {
            List<double> errorPrecursor = new List<double>(result.precursors.Count);
            foreach (Precursor precursor in result.matchedPrecursors)
            {
                PeptideSpectrumMatch psm = precursor.OptimizedBestPsm();
                if (psm.Target)
                    errorPrecursor.Add(psm.PrecursorMzError);
            }
            double variance = Numerics.Variance(errorPrecursor);
            double stdev    = Numerics.StandardDeviation(errorPrecursor);
            result.dbOptions.ConSole.WriteLine("Computed Precursor Variance = " + variance + "          STDev = " + stdev);
            if (variance < stdev)
                variance = stdev;
            //variance = result.dbOptions.precursorMassTolerance.Value * ((2 * variance) / result.dbOptions.precursorMassTolerance.Value);

            int nbRemovedPSM = 0;
            foreach (Precursor precursor in result.precursors)
            {
                for (int i = 0; i < precursor.psms.Count; )
                {
                    if (Math.Abs(precursor.psms[i].PrecursorMzError) > variance)
                    {
                        //allPSMs[i].Query.precursor.psms_AllPossibilities.Remove(allPSMs[i]);
                        precursor.psms.RemoveAt(i);
                        //allPSMs.RemoveAt(i);
                        nbRemovedPSM++;
                    }
                    else
                        i++;
                }
            }
            result.dbOptions.ConSole.WriteLine("Removed " + nbRemovedPSM + " [" + allPSMs.Count + " remaining] Peptide Spectrum matches outside the variance [" + variance + "]");
            return variance;
        }
    }
}
