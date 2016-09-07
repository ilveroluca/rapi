
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
}

// vim: set et sw=2
