
import it.crs4.rapi.*;
import it.crs4.rapi.RapiUtils;

import java.io.File;

import org.junit.*;
import static org.junit.Assert.*;

public class TestRapi
{
  private Opts rapiOpts;

  @BeforeClass
  public static void initSharedObj()
  {
    RapiUtils.loadPlugin();
  }

  @Before
  public void init() throws RapiException
  {
    rapiOpts = new Opts();
    Rapi.init(rapiOpts);
  }

  @After
  public void tearDown() throws RapiException
  {
    Rapi.shutdown();
  }

  @Test
  public void testVersions()
  {
    assertTrue(Rapi.getAlignerName().startsWith("bwa-mem"));
    assertTrue(Rapi.getAlignerVersion().startsWith("0."));
    assertTrue(Rapi.getPluginVersion().startsWith("0."));
  }

  public static void main(String args[])
  {
    TestUtils.testCaseMainMethod(TestRapi.class.getName(), args);
  }
}

// vim: set et sw=2
