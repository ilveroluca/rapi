
package it.crs4.rapi.lowrapi;

public class RapiOpNotSupportedException extends RapiException
{
  public RapiOpNotSupportedException(String reason) {
    super(reason);
  }
}
