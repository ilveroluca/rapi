
import it.crs4.rapi.lowrapi.*;
import it.crs4.rapi.RapiUtils;

import org.junit.*;
import static org.junit.Assert.*;

public class TestLowRapiRef
{
  private opts optsObj;
  private ref refObj;

  @BeforeClass
  public static void initSharedObj()
  {
    RapiUtils.loadPlugin(null);
  }

  @Before
  public void init() throws RapiException
  {
    optsObj = new opts();
    optsObj.setShare_ref_mem(false);
    Rapi.init(optsObj);
    refObj = new ref(TestUtils.RELATIVE_MINI_REF);
  }

  @After
  public void tearDown() throws RapiException
  {
    refObj.unload();
    Rapi.shutdown();
  }

  @Test
  public void testGetters()
  {
    assertEquals(1, refObj.getNContigs());
    assertTrue(refObj.getPath().length() > 1);
  }

  @Test
  public void testContig()
  {
    contig c = refObj.getContig(0);
    assertNotNull(c);

    assertEquals("chr1", c.getName());
    assertEquals(60000, c.getLen());
    assertNull(c.getAssembly_identifier());
    assertNull(c.getSpecies());
    assertNull(c.getSpecies());
    assertNull(c.getUri());
    assertNull(c.getMd5());
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testGetContigOOB1()
  {
    refObj.getContig(1);
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testGetContigOOB2()
  {
    refObj.getContig(-1);
  }

  @Test
  public void testDoubleUnload()
  {
    // double unload shouldn't cause problems
    refObj.unload();
    refObj.unload();
  }

  @Test
  public void testFormatSAMHeader() throws RapiException
  {
    String hdr = Rapi.format_sam_hdr(refObj);
    assertTrue(hdr.startsWith("@SQ"));
  }

  @Test(expected=RapiException.class)
  public void testFormatSAMHeaderError() throws RapiException
  {
    String hdr = Rapi.format_sam_hdr(null);
  }

  public static void main(String args[])
  {
    TestUtils.testCaseMainMethod(TestLowRapiRef.class.getName(), args);
  }
}

// vim: set et sw=2
