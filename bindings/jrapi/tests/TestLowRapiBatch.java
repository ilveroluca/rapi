
import it.crs4.rapi.*;
import it.crs4.rapi.RapiUtils;

import org.junit.*;
import static org.junit.Assert.*;

import java.util.List;

public class TestLowRapiBatch
{
  private Batch b;
  private List<String[]> someReads;

  @BeforeClass
  public static void initSharedObj()
  {
    RapiUtils.loadPlugin(null);
  }

  @Before
  public void init() throws RapiException
  {
    someReads = TestUtils.getSomeReads();

    Rapi.init(new Opts());
    b = new Batch(2);
  }

  @After
  public void tearDown() throws RapiException
  {
    Rapi.shutdown();
  }

  @Test
  public void testInit()
  {
    assertEquals(2, b.getNReadsPerFrag());
    assertEquals(0, b.getCapacity());
    assertEquals(0, b.getLength());
    assertEquals(0, b.getNFragments());
  }

  @Test
  public void testReserve() throws RapiException
  {
    assertEquals(0, b.getCapacity());
    b.reserve(4);
    assertEquals(2, b.getNReadsPerFrag());
    assertEquals(4, b.getCapacity());
    assertEquals(0, b.getLength());
    assertEquals(0, b.getNFragments());
  }

  @Test(expected=RapiOutOfMemoryError.class)
  public void testImpossibleReserve() throws RapiException
  {
    assertEquals(0, b.getCapacity());
    b.reserve(99999999999L);
  }


  @Test
  public void testAppend() throws RapiException
  {
    String[] aRead = someReads.get(0);

    assertEquals(2, b.getNReadsPerFrag());

    assertEquals(0, b.getLength());
    assertEquals(0, b.getNFragments());

    b.append(aRead[0], aRead[1], aRead[2], Rapi.QENC_SANGER);

    assertEquals(1, b.getLength());
    assertEquals(0, b.getNFragments());

    b.append(aRead[0], aRead[3], aRead[4], Rapi.QENC_SANGER);

    assertEquals(2, b.getLength());
    assertEquals(1, b.getNFragments());

    assertTrue(2 <= b.getCapacity());
  }

  @Test(expected=RapiOutOfMemoryError.class)
  public void testImpossibleAllocation() throws RapiException
  {
    b.reserve(2000000000);
  }

  @Test
  public void testGetRead() throws RapiException
  {
    String[] aRead = someReads.get(0);
    b.append(aRead[0], aRead[1], aRead[2], Rapi.QENC_SANGER);
    b.append(aRead[0], aRead[3], aRead[4], Rapi.QENC_SANGER);

    Read r_get = b.getRead(0, 0);

    assertNotNull(r_get);
    assertEquals(aRead[0], r_get.getId());
    assertEquals(aRead[1], r_get.getSeq());
    assertEquals(aRead[2], r_get.getQual());
    assertEquals(aRead[1].length(), r_get.getLength());

    r_get = b.getRead(0, 1);
    assertNotNull(r_get);
    assertEquals(aRead[0], r_get.getId());
    assertEquals(aRead[3], r_get.getSeq());
    assertEquals(aRead[4], r_get.getQual());
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testGetError1() throws RapiException
  {
    // verify out-of-bounds checking
    loadSomeReads(1);
    Read r = b.getRead(1, 0);
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testGetError2() throws RapiException
  {
    // verify out-of-bounds checking
    loadSomeReads(1);
    Read r = b.getRead(-1, 0);
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testGetError3() throws RapiException
  {
    // verify out-of-bounds checking
    loadSomeReads(1);
    Read r = b.getRead(0, 2);
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testGetError4() throws RapiException
  {
    // verify out-of-bounds checking
    loadSomeReads(1);
    Read r = b.getRead(0, -1);
  }

  @Test
  public void testFormatSAMBatch() throws RapiException
  {
    loadSomeReads(2);
    String output = Rapi.format_sam_batch(b);
    assertTrue(output.length() > 0);

    String[] lines = output.split("\n");
    assertEquals(4, lines.length);
    // compare Read ids.  Each element in someReads contains two reads, so two lines of SAM
    for (int i = 0; i < lines.length / 2; ++i) {
      assertEquals(someReads.get(i)[0], lines[2*i].split("\t")[0]);
      assertEquals(someReads.get(i)[0], lines[2*i+1].split("\t")[0]);
    }
  }

  @Test
  public void testFormatSAMBatchIndexed() throws RapiException
  {
    loadSomeReads(2);
    String output = Rapi.format_sam_batch(b, 1); // generate SAM for the 2nd fragment in the batch
    assertTrue(output.length() > 0);

    String[] lines = output.split("\n");
    assertEquals(2, lines.length);
    // compare Read ids
    for (int i = 0; i < lines.length; ++i) {
      assertEquals(someReads.get(1)[0], lines[i].split("\t")[0]);
    }
  }

  private void loadSomeReads(int n_fragments) throws RapiException
  {
    if (n_fragments > someReads.size())
      throw new IllegalArgumentException("Requested n_fragments (" + n_fragments + ") is greated than number of reads available for test");

    for (int i = 0; i < n_fragments; ++i) {
      String[] fragment = someReads.get(i);
      b.append(fragment[0], fragment[1], fragment[2], Rapi.QENC_SANGER);
      b.append(fragment[0], fragment[3], fragment[4], Rapi.QENC_SANGER);
    }
  }

  public static void main(String args[])
  {
    TestUtils.testCaseMainMethod(TestLowRapiBatch.class.getName(), args);
  }
}

// vim: set et sw=2
