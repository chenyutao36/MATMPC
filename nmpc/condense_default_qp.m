function mem = condense_default_qp(mem, settings)
%CONDENSE_DEFAULT_QP Store explicit full-condensing outputs.
qp = Condensing(mem, settings);
mem.Hc = qp.Hc;
mem.gc = qp.gc;
mem.Ccx = qp.Ccx;
mem.Ccg = qp.Ccg;
mem.lcc = qp.lcc;
mem.ucc = qp.ucc;
mem.lxc = qp.lxc;
mem.uxc = qp.uxc;
end
