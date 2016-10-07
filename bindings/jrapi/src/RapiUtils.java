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



package it.crs4.rapi;

import java.io.File;

public class RapiUtils
{
	private static final String SHARED_OBJ_PATH = "/jrapi.so";

  public static void loadPlugin()
  {
    loadPlugin(null);
  }

	public static void loadPlugin(String path)
	{
		if (path == null)
		{
			try {
				cz.adamh.utils.NativeUtils.loadLibraryFromJar(SHARED_OBJ_PATH);
			} catch (java.io.IOException e) {
				throw new RuntimeException(e.getMessage(), e);
			}
		}
		else
		{
			loadPluginFromPath(path);
		}
	}

  public static void loadPluginFromPath(String soPath)
  {
		if (soPath == null)
			throw new IllegalArgumentException("path must not be null");

    String absPath = new File(soPath).getAbsolutePath();
		System.err.println("Loading shared object from " + absPath);
    System.load(absPath);
  }

}

// vim: set et sw=2
