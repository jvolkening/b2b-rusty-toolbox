fn main() {
    println!("cargo:rerun-if-changed=src/alignment.cpp");
    println!("cargo:rerun-if-changed=vendor/porechop/src/alignment.cpp");
    println!("cargo:rerun-if-changed=vendor/porechop/include");
    cc::Build::new()
        .cpp(true)
        .std("c++14")
        .define("NDEBUG", None)
        .include("vendor/porechop/include")
        .file("src/alignment.cpp")
        .file("vendor/porechop/src/alignment.cpp")
        .compile("porechop_alignment");
}
