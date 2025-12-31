public class OmniPeaksStub {
    public static void main(String[] args) throws Exception {
        String prefix = "omnipeaks_stub";
        for (int i = 0; i < args.length - 1; i++) {
            if ("-p".equals(args[i])) {
                prefix = args[i + 1];
                break;
            }
        }
        java.nio.file.Files.write(java.nio.file.Paths.get(prefix + ".peak"), "chr1\t1\t2\n".getBytes());
    }
}
