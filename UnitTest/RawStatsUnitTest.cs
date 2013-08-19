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

namespace Trinity.UnitTest
{
    public static class RawStatsUnitTest
    {
        public static bool Run()
        {
            string outputPath = @"C:\_IRIC\DATA\Test\testMhc\Stats\";
            vsCSVWriter writer = new vsCSVWriter(outputPath + "output.csv");
            writer.AddLine("File,# MS1s,# MSMS,1 Charge,2 Charge,3 Charge,4 Charge,5 Charge,6 Charge,7 Charge,8 Charge,9 Charge,10 Charge,11 Charge,12 Charge,13 Charge,14 Charge");

            DBOptions options = MhcSample.CreateOptions(outputPath);
            string[] files = new string[] { @"N:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUL29_2013\Settepeptides_300713_10uL.raw", 
                                            @"N:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUL29_2013\Settepeptides_300713_10uL_MS60_MSMS15.raw",                                            
                                            @"N:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUL29_2013\Settepeptides_300713_10uL_MS60_MSMS30.raw", 
                                            @"N:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUL29_2013\Settepeptides_300713_10uL_MS60_MSMS60.raw", 
                                            @"N:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUL29_2013\Settepeptides_300713_10uL_MS120_MSMS15.raw", 
                                            @"N:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUL29_2013\Settepeptides_300713_10uL_MS120_MSMS30.raw", 
                                            @"N:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUL29_2013\Settepeptides_300713_10uL_MS120_MSMS60.raw"};
            foreach(string file in files)
            {
                pwiz.CLI.msdata.MSDataFile msFile = new pwiz.CLI.msdata.MSDataFile(file);
                Spectra spectra = Spectra.Load(msFile, options, file);
                spectra.Sort(ProductSpectrum.AscendingPrecursorMassComparison);
                
                Dictionary<Track, Precursor> DicOfComputedTracks = new Dictionary<Track,Precursor>();
                int[] charges = new int[14];
                foreach(Track track in spectra.tracks)
                    if (!DicOfComputedTracks.ContainsKey(track))
                    {
                        DicOfComputedTracks.Add(track, null);
                        int charge = 0;
                        foreach (Precursor precursor in Queries.GetIsotopes(track, options, spectra.tracks, null))
                        {
                            if (precursor.Charge > 0)
                                charge = precursor.Charge;
                            if (!DicOfComputedTracks.ContainsKey(precursor.Track))
                                DicOfComputedTracks.Add(precursor.Track, precursor);
                        }
                        charges[charge]++;
                    }
                string line = file + "," + spectra.MS1s.Count + "," + spectra.Count;
                for (int i = 0; i < charges.Length; i++)
                    line += "," + charges[i];
                writer.AddLine(line);
            }
            writer.writeToFile();
            return true;
        }
    }
}
