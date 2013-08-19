/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.IO;
using Proteomics.Utilities;

namespace Trinity
{
    /// <summary>
    /// A contaminant precursor
    /// </summary>
    public class Contaminant
    {
        public string Name;
        public string Composition;
        public double Mz;
        public Contaminant(string name, string composition, double mz)
        {
            this.Name = name;
            this.Composition = composition;
            this.Mz = mz;
        }
    }

    /// <summary>
    /// Loads potential contaminants from a tsv file
    /// </summary>
    public static class ContaminantMasses
    {
        /// <summary>
        /// Sorts from small Mz to big Mz values
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <returns></returns>
        public static int AscendingMzComparison(Contaminant left, Contaminant right)
        {
            return left.Mz.CompareTo(right.Mz);
        }

        /// <summary>
        /// List of contaminants
        /// </summary>
        private static List<Contaminant> Contaminants;

        /// <summary>
        /// Constructor that loads "contaminents.tsv" from the configuration folder
        /// </summary>
        static ContaminantMasses()
        {
            Contaminants = new List<Contaminant>();
            using (StreamReader lines = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "Configuration", "contaminents.tsv")))
            {
                string header = lines.ReadLine();

                while (lines.Peek() != -1)
                {
                    string line = lines.ReadLine();
                    string[] fields = line.Split('\t');
                    //Name  Composition Mz
                    //Methanol  CH3OH   32.0262
                    Contaminants.Add(new Contaminant(fields[0], fields[1], double.Parse(fields[2])));
                }
            }
            Contaminants.Sort(AscendingMzComparison);
            CumulContaminants = new List<int>();
            foreach (Contaminant c in Contaminants)
                CumulContaminants.Add(0);
        }

        /// <summary>
        /// Stores the number of times a contaminant was found
        /// </summary>
        public static List<int> CumulContaminants;

        /// <summary>
        /// Outputs to Console the number of time each contaminant was seen
        /// </summary>
        public static void DisplayContaminants()
        {
            Console.WriteLine("List of contaminant found ::: ");
            for (int i = 0; i < Contaminants.Count; i++)
                Console.WriteLine(Contaminants[i].Name + " " + Contaminants[i].Composition + " = " + CumulContaminants[i]);
        }

        /*
         * TODO Fix this contaminant to filter out MS values instead of MSMS fragments!!!!
        public static GraphML_List<MSPeak> RemoveContaminantsFromMzSortedList(GraphML_List<MSPeak> peaks, MassTolerance tolerance)
        {
            MSPeak[] indexes = new MSPeak[Contaminants.Count];
            double maxMz = MassTolerance.MzTop(Contaminants[Contaminants.Count - 1].Mz, tolerance);
            int c = 0;
            for(int i = 0; i < peaks.Count; i++)
            {
                MSPeak peak = peaks[i];
                if(Contaminants[c].Mz + tolerance < peak.MZ)
                    c++;

                for(int j = c; j < indexes.Length; j++)
                {
                    if (Numerics.MzDifference(peaks[i].MZ, Contaminants[j].Mz, tolerance.Units) < tolerance.Value)
                    {
                        if (indexes[j] == null || Numerics.MzDifference(peaks[i].MZ, Contaminants[j].Mz, tolerance.Units) < Numerics.MzDifference(peaks[i].MZ, indexes[j].MZ, tolerance.Units))
                            indexes[j] = peak;
                    }
                }

                if (peak.MZ > maxMz)
                    break;
            }
            for(int i = 0; i < indexes.Length; i++)
                if (indexes[i] != null)
                {
                    peaks.Remove(indexes[i]);
                    CumulContaminants[i]++;
                }
            return peaks;
        }//*/
    }
}
