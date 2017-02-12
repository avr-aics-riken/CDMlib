# CDMlib - Cartesian Data Management library


## REVISION HISTORY

---
- 2017-2-12 Version 1.0.2
  - correct install directory of tools. lib >> bin

- 2017-2-12 Version 1.0.1
  - correct FindCDM.cmake
  - enable_BUFFER_SIZE

---
- 2017-2-9 Version 1.0.0
  - cmake version
  - examples are still working
  - Tested

  |Compiler|Serial|Tools|Ex.|MPI |Tools|Ex.|
  |:--|:--:|:--:|:--:|:--:|:--:|:--:|
  |Intel 17.0.1 ||||ok|ok|0%|
  |GNU 6.2.0    |||||||
  |fx10         |||||||

---
- 2015-12-12 Version 0.9.4
  - modify a bug int tools/fconv/src/convMx1.C

---
- 2015-11-28 Version 0.9.3
  - modify assumed install directory from /usr/local/FFV to ~/FFV

---
- 2015-11-28 Version 0.9.2
  - option --with-comp need to be essential
  - change configure.ac, INSTALL, and NEWS

---
- 2015-11-27
  - fix the issue that tools/netcdf2dfi is not included in 0.9.1

---
- 2015-11-25 Version 0.9.1
  - modify INSTALL

---
- 2015-11-19 Version 0.9.0
  - Introduction of NetCDF
  - userguide_en.pdf is still previous version

---
- 2015-06-10 Version 0.8.1
  - clean package

---
- 2015-06-09 Version 0.8.0
  - introduce BUILD_DIR to keep source directory clean
  - Change to run configure
  - Change configure.ac

---
- 2015-06-05 Version 0.7.8
  - process group

---
- 2015-03-14 Version 0.7.7
  - add mpiicc, mpiicpc

---
- 2015-03-02
  - Add comment about GlobalOrigin in class cdm_Domain
  - クラス`cdm_Domain`のメンバ変数GlobalOriginに関してコメント文を追加した。

---
- 2015-02-26 Version 0.7.6
  - Change the value of argument `L_origin` of method WriteInit in convMxM and convMxN
  - FCONVのMxM変換およびMxN変換で利用している等間隔格子版のWriteInitにおいて、引数`L_origin`に各ランクの局所領域の原点座標値を与えるように変更した。

---
- 2015-02-24
  - Modify method WriteInit for uniform cartesian
  - 等間隔格子版のWriteInitの引数に各ランクの局所領域の原点座標値を与えた場合、plot3d,vtk,avs形式で出力される格子の原点座標がずれていたので、それを修正した。

---
- 2015-02-11 Version 0.7.5
  - suppress compiler warning

---
- 2015-01-31 Version 0.7.4
  - Latest autotools version

---
- 2015-01-20
  - Merge master into branch `advance_develop_2`
  - masterのリビジョン"accept pull request #13 and modify the signature of year 2015" をブランチ`advance_develop_2`にマージした。

---
- 2015-01-20
  - Bug fix for GlobalVoxel in convMxM
  - FCONVのMxM変換におけるGlobalVoxelのバグを修正した。
  - 間引き時のGlobalVoxelを求める際、領域全体のボクセル数を間引き数で割って求めるのではなく、各分割領域で間引いたボクセル数の和として求めるようにした。

---
- 2015-01-14
  - Set FortranBinary to default OutputFileType for Plot3d in FCONV
  - FCONVでPlot3d形式のファイル出力をする際、デフォルトでのファイルタイプがFortranBinaryになるように変更した。
  - また、FCONVの出力形式のチェックについて修正を加えた。

---
- 2015-01-10 Version 0.7.3
  - accept pull request #13 and modify the signature of year 2015

---
- 2015-01-09
  - Edit buffer tuning for reading BOV file and remove comments in `cdm_Array_inline.h`
  - 以前、BOV形式のファイル読み込み部をinline関数にしたが、元に戻した。
  - バッファサイズの指定をライブラリのコンパイル時に指定するようにしたので、
inline関数にする必要がなくなったため。
  - また、ファイル`cdm_Array_inline.h`内で定義しているメソッドreadBinary,writeBinaryに記載していたコメント文を削除した。

---
- 2015-01-09
  - Edit configure options for buffer tuning
  - バッファサイズの指定をライブラリのコンパイル時に指定するように、configureオプションに関する部分を修正した(元々、バッファサイズの指定はライブラリのコンパイル時ではなく、ライブラリを用いるプログラムのコンパイル時に行うようにしていた。)

---
- 2015-01-08
  - Minor corrections in conv.C and `cdm_Define.h`
  - ファイルconv.C内のif文中にあったorを||に訂正した。
  - また、ファイル`cdm_Define.h`において、エラーコードの体裁を整えた。

---
- 2015-01-05
  - Make new branch `advance_develop_2` and merge branch `advance_develop` into it
  - masterのリビジョン"accept pull request #11 and modify tools/fconv/src/Makefile.am" より新たにブランチ`advance_develop_2`を作成し、そのブランチにブランチ`advance_develop`をマージした。

---
- 2014-12-25 Version 0.7.2
  - accept pull request #11 and modify tools/fconv/src/Makefile.am

---
- 2014-12-10
  - Add message when CellID and BCflagID were converted in FCONV
  - FCONVのMxN変換およびMx1変換において、読み込んだ各ランクのCellIDとBCflagIDが
ランク0のCellIDとBCflagIDに変換される時、メッセージを出力するようにした。

---
- 2014-12-09
  - Make modifications to CellID and BCflagID
  - クラス`cdm_RankのCellID`と境界IDをランク毎に設定するように修正した。
  - また、メソッドWriteProcDfiFileにおける引数`out_host`のデフォルトをなくし、引数の順序を変更した。

---
- 2014-12-09
  - Add CellID and BCflagID in class `cdm_Rank`
  - クラス`cdm_Rank`にCellIDと境界IDを示すメンバ変数`c_id`,`bc_id`を追加し、proc.dfiのProcessのところにこれらのIDを表示するようにした。

---
- 2014-12-02 Version 0.7.1
  - add examples

---
- 2014-12-01 Version 0.7.0
  - accept pull request #10

---
- 2014-11-29 Version 0.6.3
  - Merge master into branch advance_develop
  - masterのリビジョン"accept pull request #9"をブランチadvance_developにマージした。
  - Plot3d対応で追加修正した部分を不等間隔格子対応に取り込むため。
  - その際、競合が発生したファイルが存在したので解決した。

---
- 2014-11-28
  - Minor correction about print statements in class `cdm_NonUniformDomain`
  - クラス`cdm_NonUniformDomain`内にあるprint文の訂正を行った。

---
- 2014-11-21
  - Add CoordinateFileEndian and edit class `cdm_NonUniformDomain`
  - 不等間隔格子の場合に座標ファイルのエンディアンを示す変数CoordinateFileEndianを追加した。
  - また、クラス`cdm_NonUniformDomain`のコードおよびコメント文の整理を行った。

---
- 2014-11-19 Version 0.6.2

---
- 2014-11-18
  - Edit path setting for coordinate file and add check for coordinates
  - 座標ファイルのパス指定をDFIのディレクトリからの相対パス、もしくは絶対パスで指定できるようにした。
  - また、WriteInit時に引数に与える座標配列がアロケートされているかのチェックや、ReadInit時に座標ファイル内の座標データ数のチェックを加えた。

---
- 2014-11-17
  - Merge master into branch advance_develop
  - masterのリビジョン"Add ResultFormat at VisIt section, release 0.6.1"をブランチadvance_developにマージした。
  - Plot3d対応で追加修正した部分を不等間隔格子対応に取り込むため。
  - その際、競合が発生したファイルが存在したので解決した。
  - また、Plot3d対応の追加修正において、iblankの設定方法を変更したが、単にマージしただけでは変更が反映されない部分があったので、iblankの設定方法の変更に沿って修正を加えた。

---
- 2014-11-14 Version 0.6.1

---
- 2014-11-13 Version 0.6.0

---
- 2014-11-10
  - Add check of file format for non-uniform cartesian
  - SPH形式とBOV形式は不等間隔格子に対応していないので、これらの形式で不等間隔格子の入出力を行おうとした場合はエラーを返すようにした。
  - また、不等間隔格子でproc.dfiファイルを出力する際、座標ファイルの拡張子がcrdになっているかのチェックも加えた。

---
- 2014-11-05 Version 0.5.3

---
- 2014-11-04
  - Add check for reading coordinate file
  - 不等間隔格子で座標ファイルを読み込む際、座標ファイルが読み込めているかどうか、また拡張子はcrdになっているかのチェックを付け加えた。

---
- 2014-11-02 Version 0.5.0

---
- 2014-10-30
  - Apply function overloading to method WriteInit and edit class `cdm_NonUniformDomain`
  - API関数WriteInitにオーバーロード機能の適用し、等間隔格子用と不等間隔格子用の
API関数WriteInitを用意した。
  - 等間隔格子・不等間隔格子の共通処理を行う部分は、さらにオーバーロード機能を適用したWriteInitを作成した。
  - この変更に伴い、クラスcdm_NonUniformDomainのコンストラクタ部を編集した。
  - 不等間隔格子の場合の座標データの精度は、クラスcdm_NonUniformDomainのテンプレートの型で指定するようにした。

---
- 2014-10-26
  - Minor corrections about coordinatefile information
  - 座標ファイルに関する情報の変数名やコメント文の訂正を行った。

---
- 2014-10-25 Version 0.4.0
  - Change CoordinateFileFormat to CoordinateFileType
  - 変数CoordinateFileFormatをCoordinateFileTypeに変更した。

---
- 2014-10-25
  - Merge branch advance_develop_NonUniform into branch advance_develop
  - `cdm_DFI.h`,`cdm_DFI.C`は競合が発生したので解決した。
  - `E_CDM_OUTPUT_TYPE`-> `E_CDM_FILE_TYPE`, CompName -> VariNameに変更できていない部分があったので変更した。

---
- 2014-10-24 Version 0.3.1
  - plot3d check
  - #define `D_CDM_EXT_FUNC` "func" >> "fun" // 拡張子は"fun", FileFormat識別子は"plot3d"

---
- 2014-10-22
  - Add check of number of variables for SPH format
  - SPH形式のファイルを入出力する際、変数の個数が１か３以外の場合にエラーとなるように設定した。

---
- 2014-10-22
  - Replace Component with NumVariables and edit codes following this replacement
  - dfiファイルのFileInfoの項目ComponentをNumVariablesに変更した。この変更に伴い、コード内のComponentに関する単語についてComponent->NumVariables、nComp->nVari等の置換を行った(ただし、BOVのヘッダー部のDATA_COMPONENTはそのままにした。)
  - コード内のコメントで「成分」を含んだ言葉を「変数」で置き換えた。

---
- 2014-10-22
  - Delete ArrayShape from dfi file and set ArrayShape by file format
  - dfiファイルのFileInfoの項目からArrayShapeを削除し、ArrayShapeはファイル形式によって定めるようにした。

---
- 2014-10-21
  - Restore VectorMinMax only for SPH format
  - SPH形式のみ、VectorMinMaxを入出力できるように戻した。

---
- 2014-10-21
  - Delete VectorMinMax and add check of number of variable names
  - VectorMinMaxをdfiファイルから削除した。
  - また、フィールドデータの変数の個数と登録された変数名の個数が一致するか確認するようにした。

---
- 2014-10-20
  - Set filetype in Plot3d format Fbinary
  - Plot3d形式のファイル入出力のタイプをFbinaryに設定した。

---
- 2014-10-20
  - Modify method WriteInit in class `cdm_DFI` and Read in class `cdm_FileInfo` to return error when arrayshape in plot3d is nijk
  - Plot3d形式において配列形式がnijkと指定された場合、エラーを返すようにクラス`cdm_DFI`のWriteInitとクラス`cdm_FileInfo`のReadを修正した。

---
- 2014-10-18 Version 0.2.0
  - beta version release

---
- 2014-10-17
  - Restore SPH format
  - メンバ関数ReadInitのインスタンス生成部を編集して、sph形式も読込めるように戻した。

---
- 2014-10-01
  - Change the name of `E_CDM_OUTPUT_TYPE` to `E_CDM_FILE_TYPE`
  - ファイルタイプを示す列挙型`E_CDM_OUTPUT_TYPE`の名前を`E_CDM_FILE_TYPE`に変更した。
  - クラス`cdm_DFI`のメンバ変数に`m_input_type`を追加し、この変数を用いて
plot3dで読込むファイルのタイプを指定するようにした。

---
- 2014-10-01
  - Modify methods in class `cdm_DFI_PLOT3D` and Edit method ReadInit
  - クラス`cdm_DFI_PLOT3D`内のメンバ関数の修正を行った。
  - 主に読込み部分で、バグ修正やコードの整理など。
  - 読込みに対応するファイル形式をbov,sph->bov,plot3dに変更するように、メンバ関数ReadInitのインスタンス生成部を編集した。

---
- 2014-09-26
  - Edit method `write_XYZ` in class `cdm_DFI_PLOT3D` to add iblank information to xyz file
  - Plot3dのxyzファイルにiblankの情報を追加できるようにメンバ関数`write_XYZ`を編集した。
  - iblankの情報は、メンバ関数WriteInitの引数に追加したポインタiblankを介してライブラリに取り込むようにしており、そのポインタはクラス`cdm_Domain`が保持するようにしている。

---
- 2014-09-24
  - Edit methods in class `cdm_DFI_PLOT3D` to write data including those at guide cells
  - Plot3d形式において、ガイドセル上の座標および物理量を出力するようにメンバ関数を編集した。(読込みに関しては、この時点での実装で対応できている。)

---
- 2014-09-22
  - Modify method ReadFieldData to read node data as cell-centered data
  - 格子点上のデータをセル中心のデータと見立てて読み込むように関数ReadFieldDataを修正した。

---
- 2014-09-22
  - Correct the position of the origin in xyz file
  - Plot3d形式で出力されるxyzファイルの原点の位置を訂正した。

---
- 2014-09-20
  - Correct the term written in the copyright statement in `cdm_Version.h` and `cdm_Version.h.in`
  - ファイル`cdm_Version.h`と`cdm_Version.h.in`に記載されているCopyrightの期間を訂正した。

---
- 2014-09-20
  - Change encode of file ChangeLog from Shift-JIS to UTF-8
  - ファイルChangeLogのエンコードをShift-JISからUTF-8に変更した。

---
- 2014-09-17 Version 0.1.0

---
- 2014-09-15
  - Edit method ReadFieldData to read data at nodes
  - 読み込んだデータが格子点上のデータかどうかを判断するフラグ`m_extend_arraysize_flag`を追加
  - 格子点上のデータを読み込む際、読み込み用データ配列をボクセル数より1大きい値で
作成するように関数ReadFieldDataを編集し、格子点上のデータを読み込めるようにした。
  - ヘッダーファイル`cdm_DFI_PLOT3D.h`の整理をした。

---
- 2014-09-08
  - Edit method `read_Func` in class `cdm_DFI_PLOT3D` to read Plot3d function files
  - クラス`cdm_DFI_PLOT3D`内のメンバ関数`read_Func`を用いてPlot3dファイル読み込めるように編集

---
- 2014-09-04
  - Editing class `cdm_DFI_PLOT3D`
  - クラス`cdm_DFI_PLOT3D`を編集中。
  - メンバ関数`read_Func`を用いてPlot3dファイルが読み込める状態(バッファ未使用)。

---
- 2014-09-03
  - add method read_func to class `cdm_DFI_PLOT3D` (test commit)
  - メンバ関数`read_Func`をクラス`cdm_DFI_PLOT3D`に追加。（テストコミット）

---
- 2014-09-02
  - change cio to cdm in filenames
  - ソースコード内およびソースコードのファイル名に含まれる文字cio(CIO)をcdm(CDM)に変更。

---
- 2014-09-01
  - change cio to cdm
  - ソースコード(.C, .h, .f90)内の文字cio(CIO)をcdm(CDM)に変更。

---
- 2014-08-23
  - Equivalent to CIOlib 1.5.8

---
- 2014-08-20
  - Initial commit
