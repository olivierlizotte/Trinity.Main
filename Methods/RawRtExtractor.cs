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
using System.IO;

namespace Trinity.Methods
{
    public static class RawRtExtractor
    {
        public static void FromFolderToCSV(string folder)
        {
            string[] files = Directory.GetFiles(folder, "*.RAW");
            foreach (string file in files)
            {
                ToCSV(file, folder + vsCSV.GetFileName_NoExtension(file) + "_RetentionTimes.csv");
            }
        }

        public static void ToCSV(string rawFileName, string csvOutFileName)
        {
            vsCSVWriter csvWriter = new vsCSVWriter(csvOutFileName);
            csvWriter.AddLine("Scan Number,Retention Time (min),Ms Level");
            
            pwiz.CLI.msdata.MSDataFile msFile = new pwiz.CLI.msdata.MSDataFile(rawFileName);
            
            int num_spectra = msFile.run.spectrumList.size();
            for (int i = 0; i < num_spectra; i++)
            {
                //Spectrum
                pwiz.CLI.msdata.Spectrum mySpec = msFile.run.spectrumList.spectrum(i, false);

                double retention_time = mySpec.scanList.scans[0].cvParam(pwiz.CLI.cv.CVID.MS_scan_start_time).timeInSeconds() / 60.0;
                csvWriter.AddLine((i + 1) + "," + retention_time + "," + mySpec.cvParam(pwiz.CLI.cv.CVID.MS_ms_level).value);
            }
            csvWriter.WriteToFile();
        }
    }
}
