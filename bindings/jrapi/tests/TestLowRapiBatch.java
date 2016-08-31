
import it.crs4.rapi.lowrapi.*;
import it.crs4.rapi.RapiUtils;

import org.junit.*;
import static org.junit.Assert.*;

public class TestLowRapiBatch
{
  private int error;
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
  public void init()
  {
    error = Rapi.init(new opts());
    assertEquals(Rapi.NO_ERROR, error);

    b = new batch();
    error = Rapi.reads_alloc(b, 2, 1);
  }

  @After
  public void tearDown()
  {
    Rapi.reads_free(b);
    error = Rapi.shutdown();
    assertEquals(Rapi.NO_ERROR, error);

  }

  @Test
  public void testAllocateBatch()
  {
      // batch already allocated in init()
      assertEquals(2, b.getN_reads_frag());
      assertEquals(1, b.getN_frags());
  }

  @Test
  public void testSetAndGetRead()
  {
      error = loadSomeReads(1);

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
      error = Rapi.set_read(b, 1, 0, some_reads[2][0], some_reads[2][1], some_reads[2][2], Rapi.QENC_SANGER);
      assertEquals(Rapi.NO_ERROR, error);
      error = Rapi.set_read(b, 1, 1, some_reads[3][0], some_reads[3][1], some_reads[3][2], Rapi.QENC_SANGER);
      assertEquals(Rapi.NO_ERROR, error);

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

  @Test
  public void testSetErrors()
  {
    // verify out-of-bounds checking
    error = Rapi.set_read(b, 1, 0, some_reads[0][0], some_reads[0][1], some_reads[0][2], Rapi.QENC_SANGER);
    assertFalse(error == Rapi.NO_ERROR);
    error = Rapi.set_read(b, 0, 2, some_reads[0][0], some_reads[0][1], some_reads[0][2], Rapi.QENC_SANGER);
    assertFalse(error == Rapi.NO_ERROR);
    error = Rapi.set_read(b, -1, 0, some_reads[0][0], some_reads[0][1], some_reads[0][2], Rapi.QENC_SANGER);
    assertFalse(error == Rapi.NO_ERROR);
    error = Rapi.set_read(b, 0, -1, some_reads[0][0], some_reads[0][1], some_reads[0][2], Rapi.QENC_SANGER);
    assertFalse(error == Rapi.NO_ERROR);
  }

  @Test
  public void testGetErrors()
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

  private int loadSomeReads(int n_fragments)
  {
    int error = Rapi.reads_reserve(b, n_fragments);
    if (error != Rapi.NO_ERROR)
      return error;

    int z = 0;
    for (int i = 0; i < n_fragments; ++i) {
      for (int j = 0; j < 2; ++j) {
        z = i*some_reads[0].length + j;
        error = Rapi.set_read(b, i, j, some_reads[z][0], some_reads[z][1], some_reads[z][2], Rapi.QENC_SANGER);
        if (error != Rapi.NO_ERROR)
          return error;
      }
    }
    return error;
  }

  @Test
  public void testReserve()
  {
    assertEquals(1, b.getN_frags());
    Rapi.reads_reserve(b, 1);
    assertEquals(1, b.getN_frags());
    Rapi.reads_reserve(b, 5);
    assertEquals(5, b.getN_frags());
    assertEquals(5*2, Rapi.batch_read_capacity(b));
  }
}
