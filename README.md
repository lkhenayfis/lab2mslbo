# lab2mslbo

Funcionalidades para conversão dos inputs do [`SDDPlab`](https://github.com/rjmalves/sddp-lab) para pacote [`DualSDDP`](https://github.com/bfpc/DualSDDP.jl) e testes com extensões para uso de aproximações interiores (`inner_dp`) ao invés de cortes com [`SDDP.jl`](https://github.com/odow/SDDP.jl).

## Instalação

Este repositório consiste de um projeto em `julia`, não um pacote completo, de forma que inclui o arquivo `Manifest.toml` para reprodução exata do ambiente de trabalho.

Para restauração das dependências basta realizar

```
export JULIA_SSH_NO_VERIFY_HOSTS="github.com"
export JULIA_PKG_USE_CLI_GIT=true
julia
]
activate .
instantiate
```

Caso seja necessário atualizar alguma das dependências diretas de repositórios, basta fazer, no contexto do `Pkg`de `julia`, por exemplo:

```
add https://github.com/rjmalves/sddp-lab.git#branch-tag
```

Os repositórios não listados e atualmente utilizados são o `SDDPlab` e  `DualSDDP`, que podem instalados nos branches necessários através de:

```
add https://github.com/rjmalves/sddp-lab.git#io-refactor
add https://github.com/bfpc/DualSDDP.jl.git#main
```

## Geração de imagem pré-compilada (opcional)

Algo que não é obrigatório, porém colabora para a execução de diversos testes em sequência, é a construção de imagens pré-compiladas usando o `PackageCompiler`.

`IMPORTANTE`: Somente vale realizar esta etapa caso as execuções não envolvam mudanças no código da aplicação, e sim variações dos dados de entrada. Mudanças feitas no código não serão automaticamente recompiladas se a chamada do interpretador de `julia` apontar para uma imagem pré-compilada.

```
julia
] add PackageCompiler
using PackageCompiler
] activate .
create_sysimage(; sysimage_path="lab2mslbo.so", precompile_execution_file="precompile.jl")
[ Info: PackageCompiler: Executing /home/rogerio/git/lab2mslbo/precompile.jl => /tmp/jl_packagecompiler_fgKq8t/jl_G8EEZD
[ Info: Compiling model
[ Info: Evaluating policy
...
✔ [02m:11s] PackageCompiler: compiling incremental system image
```

Para realizar a chamada do interpretador utilizando a imagem feita:

```
julia -Jlab2mslbo.so precompile.jl
```

## Uso

### SDDPlab através do DualSDDP (cálculo de política)

Um exemplo completo está disponível no arquivo `exemplo_dual.jl`, existente na raiz do repositório.

A conversão de um deck do SDDPlab se realiza através da função `build_mslbo`

```julia
include("src/lab2mslbo.jl")
M, data = lab2mslbo.build_mslbo("data-1dtoy/");
```

`M` é um objeto `DualSDDP.MSLBO` que pode ser utilizado com as demais funcionalidades do pacote e `data` é um objeto `LabData` com diversas informações necessárias para compatibilização de resultados, que foram lidas a partir do conjunto de arquivos `.jsonc` que constituem o deck do `SDDPlab`.

A partir dos objetos retornados pelo `builder`, todas as informações necessárias para se fazer uso das funcionalidades do `DualSDDP` está disponível, como as funções `primalsolve`, `primalub`, `dualsolve` e `problem_child_solve`. Além disso, foram implementadas funções adicionais para converter os arquivos de saída para um formato compatível com o utilizado pelo `SDDPlab`, facilitando a análise de resultados.

### SDDP com aproximação interior (simulação de política)

Um exemplo completo está disponível no arquivo `exemplo_inner.jl`, existente na raiz do repositório.

Para realizar a simulação considerando a aproximação interior da política (vértices ao invés de cortes), é necessário construir um objeto `SDDP.PolicyGraph` completo, contendo os vértices que representam os estados visitados.

Isto pode ser feito, por exemplo, realizando uma execução de cálculo de política normalmente, e chamando a função de conversão para um `SDDP.PolicyGraph` que possui uma `InnerBellmanFunction`. Este pode ser utilizado normalmente para a realização de simulações com `SDDP.simulate()`. No arquivo de exemplo isto é feito por meio de acesso direto a funções do `SDDPlab`, nos seus diferentes submódulos

```julia
using SDDPlab: SDDPlab
using DataFrames
using SDDP
using lab2mslbo: lab2mslbo
e = CompositeException()
entrypoint = SDDPlab.Inputs.Entrypoint("main.jsonc", e)
policy_index = findfirst(x -> isa(x, SDDPlab.Tasks.PolicyArtifact), artifacts)
policy = artifacts[policy_index].policy
inner_policy, upper_bound, upper_bound_time = lab2mslbo.__build_and_compute_ub_model(
    entrypoint.inputs.files, policy
)
```

### Cálculo de política com DualSDDP + simulação com aproximação interior

Um exemplo completo está disponível no arquivo `exemplo_dual_inner.jl`, existente na raiz do repositório.

Usando a função `read_vertices_from_file()` é possível ler um arquivo no formato `.json` exportado pelo `DualSDDP` com a política calculada, no formato de vértices para formar a aproximação interior, da maneira semelhante à do exemplo de uso anterior.

Em geral, este caso consiste em realizar uma execução que combina os dois exemplos anteriores, com uma execução do `SDDPlab` através do `DualSDDP` para cálculo de alguma política com `dualsolve()` ou `problem_child_solve()`, exportar a política para o arquivo `.json` e realizar a leitura com `read_vertices_from_file()`. O `SDDP.PolicyGraph` construído desta forma pode ser utilizado normalmente para a realização de simulações com `SDDP.simulate()`. 
