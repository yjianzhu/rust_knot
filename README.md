# rust_knot

Rust 实现的纽结拓扑分析库 — 通过 Alexander 多项式识别纽结类型并定位最小纽结核心。

本项目是对 C++ 版本（`knottype.cpp`, `knotsize.cpp`, `hull.cpp`, `my_function.cpp`）的完整重写，消除了 GiNaC 符号数学依赖，修复了已知 bug，并提供线程安全的纯函数式设计。

## 算法流程

```
输入点序列 (3D 坐标)
    │
    ▼
KMT 简化 (去除不影响拓扑的冗余点)
    │
    ▼
凸包端点延伸 (开链: 将首尾推离纽结区域)
    │
    ▼
O(n²) 交叉点检测 (XY 投影 + Z 序比较)
    │
    ▼
Alexander 矩阵构造 (Z[t] 系数)
    │
    ▼
Bareiss 无分数行列式 (替代 GiNaC)
    │
    ▼
多项式查表匹配
    │
    ▼
纽结类型 ("3_1", "4_1", ...)
    │
    ▼
二分搜索最小纽结核心 (knotsize)
    │
    ▼
KnotCoreResult { left, right, size }
```

## 项目结构

```
rust_knot/
├── Cargo.toml
├── src/
│   ├── lib.rs               # crate 入口，re-export 公开 API
│   ├── point.rs              # Point3 = [f64; 3], EPSILON = 1e-7
│   ├── error.rs              # KnotError 枚举 (Io, PolynomialParse, DataParse, HullFailed, NotFound, EmptyChain)
│   ├── polynomial.rs         # Polynomial<i64> 算术 + Bareiss 行列式 (~600 行，替代 GiNaC)
│   ├── geometry.rs           # 交叉检测、三角形相交测试、法线计算、坐标轴重定向
│   ├── alexander_table.rs    # 解析 table_knot_Alexander_polynomial.txt，HashMap 正负号双向查找
│   ├── kmt.rs                # KMT 简化算法 (开链 + 环链)
│   ├── hull.rs               # 凸包端点延伸 (基于 chull crate)
│   ├── knottype.rs           # 纽结类型识别核心: 交叉点 → Alexander 矩阵 → 多项式 → 查表
│   ├── knotsize.rs           # 二分搜索纽结核心定位
│   └── io.rs                 # XYZ / LAMMPS 格式读写
└── examples/
    └── identify_knot.rs      # CLI 工具: 识别纽结类型 + 定位核心
```

## 依赖

| crate | 版本 | 用途 |
|-------|------|------|
| `chull` | 0.2 | 3D 凸包 (QuickHull 算法) |
| `thiserror` | 2 | 错误类型派生 |
| `rayon` | 1.10 | 并行批处理 (可选使用) |
| `approx` | 0.5 | 测试用浮点比较 (dev-dependency) |

**关键设计决策**: 不使用任何符号数学库。Alexander 矩阵的元素只有 `{0, 1, -1, t, 1-t}`，因此我们实现了轻量的 `Polynomial { coeffs: Vec<i64> }`（约 600 行），配合 Bareiss 无分数消元行列式算法，完全在 `Z[t]` 整数多项式环内完成计算。

## 编译和测试

```bash
cd rust_knot

# 编译
cargo build --release

# 运行全部 36 个单元测试
cargo test

# 编译 CLI 示例
cargo build --release --example identify_knot
```

## 使用

### 作为库

```rust
use rust_knot::{AlexanderTable, get_knottype, find_knot_core, KnotOptions, KnotSizeOptions};

// 加载 Alexander 多项式查找表
let table = AlexanderTable::from_file("table_knot_Alexander_polynomial.txt")?;

// 识别纽结类型
let opts = KnotOptions { faster: true, debug: false };
let knot_type = get_knottype(&points, &table, &opts)?;  // e.g. "3_1"

// 定位最小纽结核心
let size_opts = KnotSizeOptions::default();
let core = find_knot_core(&points, &knot_type, &table, &size_opts)?;
println!("核心区间: [{}, {}], 大小: {}", core.left, core.right, core.size);
```

### CLI 工具

```bash
# 基本用法: 自动识别纽结类型并定位核心
./target/release/examples/identify_knot table_knot_Alexander_polynomial.txt input.xyz

# 指定目标纽结类型
./target/release/examples/identify_knot table_knot_Alexander_polynomial.txt input.xyz 3_1
```

输出示例:
```
Loaded Alexander table (106 entries) in 1.2ms
Read 200 points in 0.3ms
Knot type: 3_1 (computed in 5.4ms)
Knot core: [45, 162], size = 118 (found in 82.3ms)
```

## 相对 C++ 版本的改进

### Bug 修复

| Bug | C++ 位置 | Rust 修复 |
|-----|----------|-----------|
| 字符比较判断纽结复杂度 | `knotsize.cpp:259,279` — `temp[0] > target[0]` | `parse_knot_name("10_1") → (10,1)` 数值比较 |
| 未知多项式返回空字符串 | `knottype.cpp:281` | 返回 `Err(KnotError::NotFound)` |
| 全局 GiNaC mutex 阻塞并行 | `knottype.cpp:193` | 无共享可变状态，多项式运算纯函数式 |
| 凸包失败静默忽略 | `hull.cpp catch(...)` | 返回 `Option`，调用者处理回退 |
| 头文件 `using namespace std` | `myfunction.h:17` | Rust 模块系统天然隔离 |

### 架构改进

- **无 GiNaC 依赖**: 整数多项式 + Bareiss 行列式完全替代符号计算库
- **线程安全**: 所有核心函数为纯函数，无全局状态，可直接用 `rayon` 并行
- **强类型错误处理**: `KnotError` 枚举覆盖所有错误路径，无 `panic` 无静默失败
- **2530 行 Rust** 替代 ~4300 行 C++（含第三方 quickhull.h）

## 测试覆盖

36 个单元测试覆盖所有核心模块:

| 模块 | 测试数 | 覆盖内容 |
|------|--------|----------|
| `polynomial` | 16 | 四则运算、行列式、解析、归一化、精确除法、相等/哈希 |
| `geometry` | 4 | 交叉检测、法线计算、叉积 |
| `hull` | 3 | 端点延伸、退化几何、短链 |
| `kmt` | 2 | 直线简化、端点保留 |
| `knotsize` | 3 | 纽结名解析、复杂度比较、点序列旋转 |
| `io` | 4 | XYZ 读写、多帧读取、往返一致性 |
| `alexander_table` | 2 | 表解析、正负号查找 |

## 多项式查找表

程序需要 `table_knot_Alexander_polynomial.txt` 文件（位于上级目录），包含 50+ 种纽结的 Alexander 多项式。格式:

```
1                               1
3_1                             -1+t-t^2
4_1                             -1+3*t-t^2
5_1                             1-t+t^2-t^3+t^4
...
```

## License

与上级项目保持一致。
