# Como publicar pelo Git Bash

1. Extraia a pasta `arecaceae-chemical-research-landscape` dentro de:

   `G:\Meu Drive\ARTIGOS\2026\Arecaceae\github_release\`

2. Abra essa pasta, clique com o botão direito em uma área vazia e escolha **Open Git Bash here**.

3. Execute:

   ```bash
   bash PUBLISH_TO_GITHUB.sh
   ```

O script copia da pasta organizada as 15 entradas exatas, o script final, o `sessionInfo` da execução de 14/08/2026, os quatro manifests e somente os diagnósticos essenciais da etapa 05B. Em seguida, confere os checksums, exige a mesma versão do R, cria `renv.lock`, produz o manifest SHA-256, cria o commit e a tag imutável `v1.0.0` e envia tudo para a conta `carloscarollo-UFMS`.

Se o GitHub CLI (`gh`) não estiver instalado, o script fará todo o preparo local e pedirá apenas que você crie no site um repositório **público e vazio**, sem README, com este nome:

`arecaceae-chemical-research-landscape`

Depois disso, execute novamente o mesmo comando. Ele continuará do ponto em que parou.

Se a pasta do projeto estiver em outro local, informe o caminho como argumento:

```bash
bash PUBLISH_TO_GITHUB.sh "/g/Meu Drive/ARTIGOS/2026/Arecaceae"
```
