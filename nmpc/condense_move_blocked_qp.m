function mem = condense_move_blocked_qp(mem, settings)
%CONDENSE_MOVE_BLOCKED_QP Store explicit move-blocking outputs.
qp = Condensing_mb(mem, settings);
mem.Hc_r = qp.Hc_r;
mem.gc_r = qp.gc_r;
mem.Ccx_r = qp.Ccx_r;
mem.Ccg_r = qp.Ccg_r;
mem.lcc = qp.lcc;
mem.ucc = qp.ucc;
mem.lxc = qp.lxc;
mem.uxc = qp.uxc;
end
