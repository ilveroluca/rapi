
package it.crs4.rapi;

public class RapiOpNotSupportedException extends RapiException
{
  public RapiOpNotSupportedException(String reason) {
    super(reason);
  }
}
