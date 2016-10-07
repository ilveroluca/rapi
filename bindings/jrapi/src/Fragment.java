
package it.crs4.rapi;

import java.util.ArrayList;

public class Fragment implements Iterable<Read> {

  private ArrayList<Read> reads;

  public Fragment(int size) {
    reads = new ArrayList<Read>(size);
    for (int i = 0; i < size; ++i)
      reads.add(null);
  }

  public Read get(int i) {
    return reads.get(i);
  }

  public Read set(int i, Read r)
  {
    if (r == null)
      throw new IllegalArgumentException("Read cannot be null");

    return reads.set(i, r);
  }

  public java.util.Iterator<Read> iterator() {
    return reads.iterator();
  }

  public int getLength() {
    return reads.size();
  }
}
