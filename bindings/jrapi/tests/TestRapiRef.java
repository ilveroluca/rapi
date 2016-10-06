
import it.crs4.rapi.*;
import it.crs4.rapi.RapiUtils;

import org.junit.*;
import static org.junit.Assert.*;

public class TestRapiRef
{
  private Opts optsObj;
  private Ref refObj;

  @BeforeClass
  public static void initSharedObj()
  {
    RapiUtils.loadPlugin();
  }

  @Before
  public void init() throws RapiException
  {
    optsObj = new Opts();
    optsObj.setShareRefMem(false);
    Rapi.init(optsObj);
    refObj = new Ref(TestUtils.RELATIVE_MINI_REF);
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
    Contig c = refObj.getContig(0);
    assertNotNull(c);

    assertEquals("chr1", c.getName());
    assertEquals(60000, c.getLen());
    assertNull(c.getAssemblyIdentifier());
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
    String hdr = Rapi.formatSamHdr(refObj);
    assertTrue(hdr.startsWith("@SQ"));
  }

  @Test(expected=RapiException.class)
  public void testFormatSAMHeaderError() throws RapiException
  {
    String hdr = Rapi.formatSamHdr(null);
  }

  public static void main(String args[])
  {
    TestUtils.testCaseMainMethod(TestRapiRef.class.getName(), args);
  }
}

// vim: set et sw=2
