
import it.crs4.rapi.lowrapi.*;
import it.crs4.rapi.RapiUtils;

import org.junit.*;
import static org.junit.Assert.*;

public class TestLowRapiBatch
{
  private batch b;

  private static final String[][] some_reads = {
    { "read_id/1-1", "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT", "11##############################" },
    { "read_id/1-2", "TCGATCGATCGATCGATCGATCGATCGATCGA", "12$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" },
    { "read_id/2-1", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA", "21##############################" },
    { "read_id/2-2", "CGATCGATCGATCGATCGATCGATCGATCGAT", "22$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" }
  };

  @BeforeClass
  public static void initSharedObj()
  {
    RapiUtils.loadPlugin(null);
  }

  @Before
  public void init() throws RapiException
  {
    Rapi.init(new opts());

    b = new batch();
    Rapi.reads_alloc(b, 2, 1);
  }

  @After
  public void tearDown() throws RapiException
  {
    Rapi.reads_free(b);
    Rapi.shutdown();

  }

  @Test
  public void testAllocateBatch()
  {
    // batch already allocated in init()
    assertEquals(2, b.getN_reads_frag());
    assertEquals(1, b.getN_frags());
  }

  @Test(expected=RapiOutOfMemoryError.class)
  public void testImpossibleAllocation() throws RapiException
  {
    b = new batch();
    Rapi.reads_alloc(b, 2, 2000000000);
  }

  @Test
  public void testSetAndGetRead() throws RapiException
  {
    loadSomeReads(1);

    read r_get = Rapi.get_read(b, 0, 0);
    assertNotNull(r_get);
    assertEquals(some_reads[0][0], r_get.getId());
    assertEquals(some_reads[0][1], r_get.getSeq());
    assertEquals(some_reads[0][2], r_get.getQual());
    assertEquals(some_reads[0][1].length(), r_get.getLength());

    r_get = Rapi.get_read(b, 0, 1);
    assertNotNull(r_get);
    assertEquals(some_reads[1][0], r_get.getId());
    assertEquals(some_reads[1][1], r_get.getSeq());
    assertEquals(some_reads[1][2], r_get.getQual());

    Rapi.reads_reserve(b, 2);
    Rapi.set_read(b, 1, 0, some_reads[2][0], some_reads[2][1], some_reads[2][2], Rapi.QENC_SANGER);
    Rapi.set_read(b, 1, 1, some_reads[3][0], some_reads[3][1], some_reads[3][2], Rapi.QENC_SANGER);

    r_get = Rapi.get_read(b, 1, 0);
    assertNotNull(r_get);
    assertEquals(some_reads[2][0], r_get.getId());
    assertEquals(some_reads[2][1], r_get.getSeq());
    assertEquals(some_reads[2][2], r_get.getQual());

    r_get = Rapi.get_read(b, 1, 1);
    assertNotNull(r_get);
    assertEquals(some_reads[3][0], r_get.getId());
    assertEquals(some_reads[3][1], r_get.getSeq());
    assertEquals(some_reads[3][2], r_get.getQual());
  }

  // verify out-of-bounds checking
  @Test(expected=RapiInvalidParamException.class)
  public void testSetErrors1() throws RapiException {
    Rapi.set_read(b, 1, 0, some_reads[0][0], some_reads[0][1], some_reads[0][2], Rapi.QENC_SANGER);
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testSetErrors2() throws RapiException {
    Rapi.set_read(b, 0, 2, some_reads[0][0], some_reads[0][1], some_reads[0][2], Rapi.QENC_SANGER);
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testSetErrors3() throws RapiException {
    Rapi.set_read(b, -1, 0, some_reads[0][0], some_reads[0][1], some_reads[0][2], Rapi.QENC_SANGER);
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testSetErrors4() throws RapiException {
    Rapi.set_read(b, 0, -1, some_reads[0][0], some_reads[0][1], some_reads[0][2], Rapi.QENC_SANGER);
  }

  @Test
  public void testGetErrors() throws RapiException
  {
    loadSomeReads(1);
    // verify out-of-bounds checking
    read r = Rapi.get_read(b, 0, 0);
    assertNotNull(r);
    r = Rapi.get_read(b, 1, 0);
    assertNull(r);
    r = Rapi.get_read(b, 0, 2);
    assertNull(r);
    r = Rapi.get_read(b, -1, 0);
    assertNull(r);
    r = Rapi.get_read(b, 0, -1);
    assertNull(r);
  }

  @Test
  public void testFormatSAMBatch() throws RapiException
  {
    loadSomeReads(2);
    String output = Rapi.format_sam_batch(b);
    assertTrue(output.length() > 0);

    String[] lines = output.split("\n");
    assertEquals(4, lines.length);
  }


  private void loadSomeReads(int n_fragments) throws RapiException
  {
    if (n_fragments > some_reads.length / 2)
      throw new IllegalArgumentException("Requested n_fragments (" + n_fragments + ") is greated than number of reads available for test");

    Rapi.reads_reserve(b, n_fragments);

    int z = 0;
    for (int i = 0; i < n_fragments; ++i) {
      for (int j = 0; j < 2; ++j) {
        z = i*2 + j;
        Rapi.set_read(b, i, j, some_reads[z][0], some_reads[z][1], some_reads[z][2], Rapi.QENC_SANGER);
      }
    }
  }

  @Test
  public void testReserve() throws RapiException
  {
    assertEquals(1, b.getN_frags());
    Rapi.reads_reserve(b, 1);
    assertEquals(1, b.getN_frags());
    Rapi.reads_reserve(b, 5);
    assertEquals(5, b.getN_frags());
    assertEquals(5*2, Rapi.batch_read_capacity(b));
  }
}

// vim: set et sw=2
