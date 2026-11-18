"""用 Windows 作业对象限制当前计算进程及其后续子进程的总资源。"""

import ctypes
import os
from ctypes import wintypes


class BasicLimits(ctypes.Structure):
    _fields_ = [
        ("PerProcessUserTimeLimit", ctypes.c_longlong), ("PerJobUserTimeLimit", ctypes.c_longlong),
        ("LimitFlags", wintypes.DWORD), ("MinimumWorkingSetSize", ctypes.c_size_t),
        ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", wintypes.DWORD),
        ("Affinity", ctypes.c_size_t), ("PriorityClass", wintypes.DWORD), ("SchedulingClass", wintypes.DWORD),
    ]


class IOCounters(ctypes.Structure):
    _fields_ = [(name, ctypes.c_ulonglong) for name in ["ReadOperationCount", "WriteOperationCount", "OtherOperationCount", "ReadTransferCount", "WriteTransferCount", "OtherTransferCount"]]


class ExtendedLimits(ctypes.Structure):
    _fields_ = [("BasicLimitInformation", BasicLimits), ("IoInfo", IOCounters)] + [(name, ctypes.c_size_t) for name in ["ProcessMemoryLimit", "JobMemoryLimit", "PeakProcessMemoryUsed", "PeakJobMemoryUsed"]]


class CPURate(ctypes.Structure):
    _fields_ = [("ControlFlags", wintypes.DWORD), ("CpuRate", wintypes.DWORD)]


def apply_limits(memory_bytes=6_000_000_000, cpu_percent=50):
    if os.name != "nt":
        raise RuntimeError("本资源限制器需要 Windows 作业对象")
    if not 0 < memory_bytes <= 8_000_000_000 or not 0 < cpu_percent <= 50:
        raise ValueError("不得超出用户最新授权的8 GB和50% CPU上限")
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
    kernel.CreateJobObjectW.restype = wintypes.HANDLE
    kernel.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
    kernel.QueryInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD, ctypes.c_void_p]
    kernel.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    kernel.CloseHandle.argtypes = [wintypes.HANDLE]
    handle = kernel.CreateJobObjectW(None, None)
    if not handle:
        raise OSError(ctypes.get_last_error(), "创建资源限制对象失败")

    def require(value, message):
        if not value:
            error = ctypes.get_last_error()
            kernel.CloseHandle(handle)
            raise OSError(error, message)

    limits = ExtendedLimits()
    limits.BasicLimitInformation.LimitFlags = 0x200
    limits.JobMemoryLimit = memory_bytes
    require(kernel.SetInformationJobObject(handle, 9, ctypes.byref(limits), ctypes.sizeof(limits)), "设置总内存硬限失败")
    cpu = CPURate(0x1 | 0x4, int(cpu_percent * 100))
    require(kernel.SetInformationJobObject(handle, 15, ctypes.byref(cpu), ctypes.sizeof(cpu)), "设置 CPU 硬限失败")
    # 先限制控制器，再启动工作进程，后者自动继承同一作业的总配额。
    require(kernel.AssignProcessToJobObject(handle, kernel.GetCurrentProcess()), "加入资源限制对象失败")

    def snapshot():
        current = ExtendedLimits()
        rate = CPURate()
        if not kernel.QueryInformationJobObject(handle, 9, ctypes.byref(current), ctypes.sizeof(current), None):
            raise OSError(ctypes.get_last_error(), "查询内存配额失败")
        if not kernel.QueryInformationJobObject(handle, 15, ctypes.byref(rate), ctypes.sizeof(rate), None):
            raise OSError(ctypes.get_last_error(), "查询 CPU 配额失败")
        return {
            "memory_hard_limit_bytes": current.JobMemoryLimit,
            "cpu_hard_limit_percent": rate.CpuRate / 100,
            "cpu_control_flags": rate.ControlFlags,
            "includes_controller_and_worker": True,
            "peak_job_committed_memory_bytes": current.PeakJobMemoryUsed,
        }

    return snapshot, lambda: kernel.CloseHandle(handle)
