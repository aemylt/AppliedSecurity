import java.util.Arrays;
import javax.crypto.Cipher;
import javax.crypto.spec.SecretKeySpec;
 
public class encrypt {
  public static void main( String[] args ) {
    try {
      byte[] k = new byte[] { ( byte )( 0x2B ), ( byte )( 0x7E ), 
                              ( byte )( 0x15 ), ( byte )( 0x16 ), 
                              ( byte )( 0x28 ), ( byte )( 0xAE ), 
                              ( byte )( 0xD2 ), ( byte )( 0xA6 ),
                              ( byte )( 0xAB ), ( byte )( 0xF7 ), 
                              ( byte )( 0x15 ), ( byte )( 0x88 ), 
                              ( byte )( 0x09 ), ( byte )( 0xCF ), 
                              ( byte )( 0x4F ), ( byte )( 0x3C ) };
      byte[] m = new byte[] { ( byte )( 0x32 ), ( byte )( 0x43 ), 
                              ( byte )( 0xF6 ), ( byte )( 0xA8 ), 
                              ( byte )( 0x88 ), ( byte )( 0x5A ), 
                              ( byte )( 0x30 ), ( byte )( 0x8D ),
                              ( byte )( 0x31 ), ( byte )( 0x31 ), 
                              ( byte )( 0x98 ), ( byte )( 0xA2 ), 
                              ( byte )( 0xE0 ), ( byte )( 0x37 ), 
                              ( byte )( 0x07 ), ( byte )( 0x34 ) };
      byte[] c = new byte[] { ( byte )( 0x39 ), ( byte )( 0x25 ), 
                              ( byte )( 0x84 ), ( byte )( 0x1D ), 
                              ( byte )( 0x02 ), ( byte )( 0xDC ), 
                              ( byte )( 0x09 ), ( byte )( 0xFB ),
                              ( byte )( 0xDC ), ( byte )( 0x11 ), 
                              ( byte )( 0x85 ), ( byte )( 0x97 ), 
                              ( byte )( 0x19 ), ( byte )( 0x6A ), 
                              ( byte )( 0x0B ), ( byte )( 0x32 ) };
      byte[] t;

      Cipher cipher = Cipher.getInstance( "AES/ECB/NoPadding", "SunJCE" );
      cipher.init( Cipher.ENCRYPT_MODE, new SecretKeySpec( k, "AES" ) );
      t = cipher.doFinal( m );

      if( Arrays.equals( t, c ) ) {
        System.out.println( "AES.Enc( k, m ) == c" );
      }
      else {
        System.out.println( "AES.Enc( k, m ) != c" );
      }
    } catch( Exception e ) {
      e.printStackTrace();
    } 
  }
}
