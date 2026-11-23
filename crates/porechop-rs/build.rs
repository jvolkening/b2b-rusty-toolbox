fn main() {
    println!("cargo:rerun-if-changed=src/alignment.cpp");
    println!("cargo:rerun-if-changed=src/middle_alignment.cpp");
    println!("cargo:rerun-if-changed=vendor/porechop/src/alignment.cpp");
    println!("cargo:rerun-if-changed=vendor/porechop/include");
    println!("cargo:rerun-if-changed=vendor/parasail");
    let target = std::env::var("CARGO_CFG_TARGET_ARCH").expect("target architecture");
    let x86 = target == "x86" || target == "x86_64";
    let include = "vendor/parasail";
    let mut native = cc::Build::new();
    native
        .include(include)
        .define("_POSIX_C_SOURCE", Some("200112L"))
        .warnings(false)
        .files([
            "vendor/parasail/src/memory.c",
            "vendor/parasail/src/sg_trace.c",
            "vendor/parasail/src/traceback.c",
        ]);
    native.compile("parasail_common");
    if x86 {
        let msvc = std::env::var("CARGO_CFG_TARGET_ENV").as_deref() == Ok("msvc");
        let mut sse = cc::Build::new();
        sse.include(include).warnings(false).files([
            "vendor/parasail/src/memory_sse.c",
            "vendor/parasail/src/sg_trace_striped_sse2_128_32.c",
        ]);
        if msvc {
            if target == "x86" {
                sse.flag("/arch:SSE2");
            }
        } else {
            sse.flag("-msse2");
        }
        sse.compile("parasail_sse2");
        let mut avx = cc::Build::new();
        avx.include(include).warnings(false).files([
            "vendor/parasail/src/memory_avx2.c",
            "vendor/parasail/src/sg_trace_striped_avx2_256_32.c",
        ]);
        if msvc {
            avx.flag("/arch:AVX2");
        } else {
            avx.flag("-mavx2");
        }
        avx.compile("parasail_avx2");
    }
    cc::Build::new()
        .cpp(true)
        .std("c++14")
        .define("NDEBUG", None)
        .include("vendor/porechop/include")
        .include(include)
        .file("src/alignment.cpp")
        .file("src/middle_alignment.cpp")
        .file("vendor/porechop/src/alignment.cpp")
        .compile("porechop_alignment");
}
