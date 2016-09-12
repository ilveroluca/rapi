
import it.crs4.rapi.lowrapi.*;
import it.crs4.rapi.RapiUtils;

import java.io.File;

import org.junit.*;
import static org.junit.Assert.*;

public class TestLowRapi
{
  private opts rapiOpts;

  @BeforeClass
  public static void initSharedObj()
  {
    RapiUtils.loadPlugin(null);
  }

  @Before
  public void init() throws RapiException
  {
    rapiOpts = new opts();
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
    assertTrue(Rapi.aligner_name().startsWith("bwa-mem"));
    assertTrue(Rapi.aligner_version().startsWith("0."));
    assertTrue(Rapi.plugin_version().startsWith("0."));
  }

  public static void main(String args[])
  {
    TestUtils.testCaseMainMethod(TestLowRapi.class.getName(), args);
  }
}

// vim: set et sw=2
