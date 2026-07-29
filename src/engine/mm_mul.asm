; Hand-written x86-64 (MASM/ml64) big-integer kernels for the Mandelbrot
; reference orbit. De-risking step: a mulx-based mpn_addmul_1 equivalent to see
; how close hand-asm gets to GMP's tuned kernel. Windows x64 ABI:
;   rcx=arg0 rdx=arg1 r8=arg2 r9=arg3 ; callee-saved: rbx rbp rsi rdi r12-r15.

OPTION CASEMAP:NONE
.code

; uint64_t mm_addmul_1(uint64_t* rp, const uint64_t* up, int64_t n, uint64_t v)
;   rp[0..n-1] += up[0..n-1]*v ; returns the carry out of the top limb.
;   Single carry chain, mulx for the multiply. n assumed >= 1.
mm_addmul_1 PROC
    push    rbx
    mov     r10, rdx            ; r10 = up
    mov     rdx, r9             ; rdx = v (mulx implicit multiplier)
    xor     r11d, r11d          ; j = 0
    xor     r9d, r9d            ; carry = 0
Laddmul_loop:
    mulx    rbx, rax, QWORD PTR [r10 + r11*8]   ; rbx:rax = v * up[j]
    add     rax, r9                             ; lo += carry
    adc     rbx, 0                              ; hi += CF
    add     QWORD PTR [rcx + r11*8], rax        ; rp[j] += lo
    adc     rbx, 0                              ; hi += CF
    mov     r9, rbx                             ; carry = hi
    inc     r11
    cmp     r11, r8
    jl      Laddmul_loop
    mov     rax, r9             ; return carry
    pop     rbx
    ret
mm_addmul_1 ENDP

END
