using Statistics
using QuantumClifford
using Random
using StableRNGs; 
rng = StableRNG(77);

function ZXZstate(N_sites);
    # ZXZ state
    p0 = zeros(N_sites);p0 = convert.(Bool,p0);
    x1 = copy(p0);x1[1]=1;z1 = copy(p0);z1[2] = 1;z1[end] = 1;
    Pind = PauliOperator(0x0,x1,z1);
    for k in 2:N_sites-1
        x1 = copy(p0);x1[k] = 1;
        z1 = copy(p0);z1[k-1] = 1;z1[k+1] = 1;
        zxz = PauliOperator(0x0,x1,z1);
        Pind = [Pind;zxz]
    end
    x1 = copy(p0);x1[end]=1;z1 = copy(p0);z1[1] = 1;z1[N_sites-1] = 1;
    Pind =[Pind;PauliOperator(0x0,x1,z1)];
    ψ_zxz = Stabilizer(Pind);
    return ψ_zxz
end

N_sites = 20;Sample = 1e5;ψ_ghz = ghz(N_sites);
ψ_zxz = ZXZstate(N_sites);
Omk = [];Ovk = [];k_sub = collect(5:15);
for ks in 1:length(k_sub)
    k_regim = k_sub[ks];
    n1 = 1;n2 = k_regim;
    Oave = [];state = copy(ψ_zxz);#state = copy(ψ_ghz);
    x1 = zeros(N_sites);x1 = convert.(Bool,x1);
    #z1 = copy(x1);z1[n1:n2] .= Bool[1];
    #Z₁...Zₖ string operator for GHZ state.
    z1 = copy(x1);z1[n1:n1+1] .= 1;z1[n2-1:n2] .= 1;x1[n1+1:n2-1] .= Bool[1]
    #Z₁Y₂X₃... string operator for ZXZ state.
    Zmeasure = PauliOperator(0x0,x1,z1);
    om = [];
    for saa in 1:Sample
        a=random_clifford(rng,k_regim);ai = inv(a);
        b=apply!(copy(state),a,collect(n1:n2));
        bmix = MixedDestabilizer(b)
        for k in n1:n2
            projectZrand!(bmix,k)
        end
        bmix = stabilizerview(bmix);
        cn = apply!(copy(bmix),ai,collect(n1:n2));
        Omeas = expect(Zmeasure,cn);
        push!(Oave,Omeas);
    end
    Wₕₐ = 2^k_regim+1;# shadow norm of random Clifford ensemble.
    Oave .*=Wₕₐ;
    Oest = mean(Oave);
    Oacc = (-1)^k_regim;
    Ovar = mean(Oave .^2)-Oest.^2;
    push!(Omk,Oest);push!(Ovk,Ovar);
end

Wei_Haar = 2 .^(k_sub) .+1;
