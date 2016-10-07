/************************************************************************************
 * This code is published under the The MIT License.
 *
 * Copyright (c) 2016 Center for Advanced Studies,
 *                      Research and Development in Sardinia (CRS4), Pula, Italy.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ************************************************************************************/



import it.crs4.rapi.RapiUtils;
import it.crs4.rapi.RapiException;
import it.crs4.rapi.RapiConstants;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TestUtils
{
  public static final String RELATIVE_MINI_REF = "../../tests/mini_ref/mini_ref.fasta";
  public static final String RELATIVE_MINI_REF_SEQS = "../../tests/mini_ref/mini_ref_seqs.txt";

  public static final String[][] some_reads = {
    { "read_id:1", "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT", "11##############################",
                     "TCGATCGATCGATCGATCGATCGATCGATCGA", "12$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" },
    { "read_id:2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA", "21##############################",
                     "CGATCGATCGATCGATCGATCGATCGATCGAT", "22$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" }
  };


  public static List<String[]> getSomeReads()
  {
    List<String[]> list = new ArrayList<String[]>(some_reads.length);
    for (int i = 0; i < some_reads.length; ++i)
      list.add(some_reads[i]);

    return list;
  }

  /**
   *  Fetch a tab-separated file of sequences::
   *
   *      ID	SEQ1	Q1	SEQ2	Q2
   */

  public static List<String[]> readSeqs(String path) throws FileNotFoundException, IOException
  {
    List<String[]> list = new ArrayList<String[]>(5);
    try(BufferedReader br = new BufferedReader(new FileReader(path)))
    {
      String line = br.readLine();
      while (line != null)
      {
        String[] record = line.split("\t");
        if (record.length != 5) {
          throw new RuntimeException("Invalid read file format.  Found " + record.length +
              " tab-separated fields while we expected 5");
        }
        list.add(record);
        line = br.readLine();
      }
    }
    return list;
  }

  public static List<String[]> readMiniRefSeqs() throws FileNotFoundException, IOException
  {
    return readSeqs(RELATIVE_MINI_REF_SEQS);
  }


  public static void appendSeqsToBatch(List<String[]> seqs, it.crs4.rapi.Batch batch) throws RapiException
  {
    if (batch.getNReadsPerFrag() != 2)
      throw new RuntimeException("This method assumes batches of paired reads (getNReadsPerFrag == 2)");

    for (String[] fragment: seqs)
    {
      batch.append(fragment[0], fragment[1], fragment[2], RapiConstants.QENC_SANGER);
      batch.append(fragment[0], fragment[3], fragment[4], RapiConstants.QENC_SANGER);
    }
  }


  public static void testCaseMainMethod(String testName, String[] args)
  {
    if (args.length == 1)
      RapiUtils.loadPlugin(args[0]);
    else
      RapiUtils.loadPlugin();

    org.junit.runner.JUnitCore.main(testName);
  }
}

// vim: set et sw=2
