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
