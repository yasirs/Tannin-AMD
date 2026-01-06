# GPU / CUDA / PyTorch Environment Summary

**Last verified:** 2026-01-06  
**Host:** tag-737  
**OS:** Ubuntu 24.04.3 LTS  
**Kernel:** 6.8.0-90-generic  

This file documents the GPU, CUDA, and PyTorch environment so that
automation tools and LLM agents can reason about available capabilities
without re-running system discovery.

---

## Hardware

- **GPUs:** 2 × NVIDIA L40 (Ada Lovelace)
- **Architecture:** AD102
- **Compute Capability:** sm_89
- **Memory:** ~48 GB per GPU
- **NVLink:** ❌ No (PCIe only; same NUMA node)
- **Power Limit:** 300 W (max, default)
- **Persistence Mode:** Enabled
- **ECC:** Supported (driver-managed)

Topology summary:
- GPU0 ↔ GPU1 connection type: `NODE` (PCIe, same NUMA)
- CPU affinity: 0–127
- NUMA node: 0

---

## NVIDIA Driver

- **Driver Version:** 580.105.08
- **Kernel Module:** Open kernel module (DKMS)
- **CUDA Driver API Support:** 13.0
- **Status:** Verified healthy via `nvidia-smi`

---

## CUDA Toolkit

- **CUDA Toolkit Version:** 12.9
- **nvcc:** 12.9.86
- **Install Method:** Apt-managed (NOT runfile)
- **CUDA Home:** `/usr/local/cuda` → `/usr/local/cuda-12.9`
- **Runtime Libraries:** Found via `ldconfig`
- **cuDNN:** cuDNN 9 (CUDA 12.x compatible)

Custom CUDA code compilation is verified to work:
- `nvcc` compilation succeeds
- Kernel launch + device sync succeeds on L40

Recommended architecture flags for custom builds:
- `-arch=sm_89`
- or `-gencode arch=compute_89,code=sm_89`

---

## Python / PyTorch

- **Python:** 3.12.3
- **Environment:** `/home/ysuhail/venvs/torch`
- **PyTorch:** `2.8.0+cu129`
- **torch.version.cuda:** 12.9
- **CUDA Available:** ✅ True
- **Device Count:** 2
- **NCCL:** Available and functional
- **cuDNN:** Available

Verified via:
- GPU tensor allocation
- fp16 matrix multiplication
- Observed GPU memory usage in `nvidia-smi`
- NCCL availability check

---

## Multi-GPU Notes

- Multi-GPU training (DDP) is supported via NCCL.
- No NVLink → expect PCIe bandwidth, not NVLink-class.
- Suitable for:
  - DataParallel / DDP
  - Model parallelism (with PCIe-aware assumptions)
  - Large-batch training and inference

---

## Known-Good State

As of last verification:
- NVIDIA driver loads correctly after reboot
- CUDA toolkit and runtime are consistent
- PyTorch CUDA build matches system CUDA
- No environment variables masking GPUs (`CUDA_VISIBLE_DEVICES` unset)

If CUDA or PyTorch stop working:
1. Check `nvidia-smi`
2. Check running kernel vs installed kernel
3. Reconfirm active Python environment
4. Verify `torch.version.cuda` vs driver

---

## Do NOT Do

- Do not install NVIDIA `.run` drivers
- Do not downgrade drivers below 550 on Ubuntu 24.04
- Do not install multiple CUDA toolkits unless explicitly required
- Do not assume NVLink is present

---

## Intended Workloads

This system is suitable for:
- PyTorch training/inference (fp16/bf16)
- Custom CUDA / C++ extensions
- Multi-GPU experiments via NCCL
- LLM workloads, scRNA-seq models, VAEs, transformers
- torch.compile / modern PyTorch features

---
