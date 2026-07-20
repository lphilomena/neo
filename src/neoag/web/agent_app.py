import json, os, re, signal, subprocess, time, uuid
from pathlib import Path
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse, PlainTextResponse
from pydantic import BaseModel

ROOT = Path(os.environ.get('NEOAG_PROJECT_ROOT', Path.cwd())).resolve()
JOBS_DIR = Path(os.environ.get('NEOAG_AGENT_WEB_JOBS', ROOT / 'results' / 'agent_web_jobs')).resolve()
UPLOADS_DIR = Path(os.environ.get('NEOAG_AGENT_WEB_UPLOADS', ROOT / 'results' / 'agent_web_uploads')).resolve()
JOBS_DIR.mkdir(parents=True, exist_ok=True)
UPLOADS_DIR.mkdir(parents=True, exist_ok=True)
app = FastAPI(title='NeoAg LLM Agent Web')
JOBS = {}

class Req(BaseModel):
    message: str
    files: list[str] = []
    outdir: str = 'results/llm_agent_web'
    mode: str = 'execute-with-approval'
    llm_provider: str = 'deepseek'
    model: str = 'deepseek-chat'
    api_key_env: str = 'DEEPSEEK_API_KEY'
    api_base: str | None = None
    sample_id: str | None = None
    allow_high_risk: bool = True

def jp(job_id):
    d = JOBS_DIR / job_id
    return d, d / 'job.json', d / 'agent.log'

def safe_name(value, default='upload'):
    name = Path(str(value or default)).name
    name = re.sub(r'[^A-Za-z0-9._+-]+', '_', name).strip('._')
    return name or default

def save(job):
    d, meta, _ = jp(job['job_id'])
    d.mkdir(parents=True, exist_ok=True)
    meta.write_text(json.dumps({k: v for k, v in job.items() if k != 'proc'}, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

def load(job_id):
    if job_id in JOBS:
        return JOBS[job_id]
    _, meta, _ = jp(job_id)
    if not meta.exists():
        raise HTTPException(404, 'job not found')
    job = json.loads(meta.read_text(encoding='utf-8'))
    JOBS[job_id] = job
    return job

def pub(job):
    proc = job.get('proc')
    if proc:
        rc = proc.poll()
        if rc is None:
            job['status'] = 'RUNNING'
        else:
            job.update(status='PASS' if rc == 0 else 'FAIL', returncode=rc, finished_at=time.time())
            job.pop('proc', None)
            save(job)
    elif job.get('status') == 'RUNNING' and job.get('pid'):
        pid_alive = Path('/proc') .joinpath(str(job['pid'])).exists()
        if not pid_alive:
            log_path = jp(job['job_id'])[2]
            log_text = log_path.read_text(encoding='utf-8', errors='ignore') if log_path.exists() else ''
            if '"status": "PASS"' in log_text or "'status': 'PASS'" in log_text:
                job.update(status='PASS', returncode=0, finished_at=job.get('finished_at') or time.time())
            elif '"status": "FAIL"' in log_text or 'Traceback' in log_text or 'ERROR' in log_text:
                job.update(status='FAIL', returncode=1, finished_at=job.get('finished_at') or time.time())
            else:
                job.update(status='STOPPED', returncode=-15, finished_at=job.get('finished_at') or time.time())
            save(job)
    out = {k: v for k, v in job.items() if k != 'proc'}
    out['log_path'] = str(jp(job['job_id'])[2])
    return out

@app.get('/', response_class=HTMLResponse)
def home():
    return HTMLResponse(HTML, headers={'Cache-Control': 'no-store, max-age=0', 'Pragma': 'no-cache'})

@app.get('/health')
def health():
    return {'ok': True, 'project_root': str(ROOT), 'jobs_dir': str(JOBS_DIR)}

@app.post('/api/upload')
async def upload(request: Request, filename: str = 'upload.dat', upload_id: str | None = None):
    group = safe_name(upload_id or (time.strftime('%Y%m%d_%H%M%S') + '_' + uuid.uuid4().hex[:8]), 'upload')
    name = safe_name(filename, 'upload.dat')
    target_dir = UPLOADS_DIR / group
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / name
    size = 0
    with target.open('wb') as fh:
        async for chunk in request.stream():
            if chunk:
                fh.write(chunk)
                size += len(chunk)
    if size == 0:
        target.unlink(missing_ok=True)
        raise HTTPException(400, 'empty upload')
    return {'path': str(target), 'size': size, 'upload_id': group}

@app.get('/api/jobs')
def jobs():
    rows = []
    for meta in sorted(JOBS_DIR.glob('*/job.json'), key=lambda p: p.stat().st_mtime, reverse=True):
        try:
            rows.append(pub(load(meta.parent.name)))
        except Exception:
            pass
    return {'jobs': rows[:100]}

@app.post('/api/jobs')
def run(req: Req):
    if req.mode not in {'plan', 'execute-safe', 'execute-with-approval'}:
        raise HTTPException(400, 'bad mode')
    job_id = time.strftime('%Y%m%d_%H%M%S') + '_' + uuid.uuid4().hex[:8]
    d, _, log = jp(job_id)
    d.mkdir(parents=True, exist_ok=True)
    cmd = [str(ROOT / 'bin' / 'neoag-llm-agent'), '--llm-provider', req.llm_provider, '--model', req.model, '--api-key-env', req.api_key_env, '--message', req.message, '--outdir', req.outdir, '--mode', req.mode, '--project-root', str(ROOT)]
    if req.allow_high_risk:
        cmd.append('--allow-high-risk')
    if req.api_base:
        cmd += ['--api-base', req.api_base]
    if req.sample_id:
        cmd += ['--sample-id', req.sample_id]
    for f in req.files:
        if f.strip():
            cmd += ['--file', f.strip()]
    env = os.environ.copy()
    env['PYTHONPATH'] = f"{ROOT / 'src'}:{env.get('PYTHONPATH', '')}"
    env['PYTHONUNBUFFERED'] = '1'
    with log.open('w', encoding='utf-8') as fh:
        proc = subprocess.Popen(cmd, cwd=str(ROOT), env=env, stdout=fh, stderr=subprocess.STDOUT, text=True, start_new_session=True)
    job = {'job_id': job_id, 'status': 'RUNNING', 'pid': proc.pid, 'created_at': time.time(), 'outdir': req.outdir, 'message': req.message, 'files': req.files, 'command': cmd, 'proc': proc}
    JOBS[job_id] = job
    save(job)
    return pub(job)

@app.get('/api/jobs/{job_id}')
def one(job_id):
    return pub(load(job_id))

@app.post('/api/jobs/{job_id}/stop')
def stop(job_id):
    job = load(job_id)
    if job.get('pid'):
        try:
            os.killpg(int(job['pid']), signal.SIGTERM)
        except Exception:
            pass
    job.update(status='STOPPED', returncode=-15, finished_at=time.time())
    job.pop('proc', None)
    save(job)
    return pub(job)

@app.get('/api/jobs/{job_id}/log', response_class=PlainTextResponse)
def get_log(job_id, tail: int = 500):
    load(job_id)
    p = jp(job_id)[2]
    if not p.exists():
        return ''
    lines = p.read_text(encoding='utf-8', errors='ignore').splitlines()
    return '\n'.join(lines[-tail:]) + ('\n' if lines else '')

@app.get('/api/jobs/{job_id}/final', response_class=PlainTextResponse)
def get_final(job_id):
    job = load(job_id)
    p = ROOT / job.get('outdir', '') / 'final_response.md'
    if not p.exists():
        return ''
    created_at = float(job.get('created_at') or 0)
    # Multiple web jobs may share the same outdir. Do not show a previous job's
    # final_response.md while the current job is still generating its answer.
    if created_at and p.stat().st_mtime < created_at:
        return ''
    return p.read_text(encoding='utf-8', errors='ignore')

@app.get('/api/jobs/{job_id}/outputs')
def outs(job_id):
    job = load(job_id)
    outdir = ROOT / job.get('outdir', '')
    rels = ['case_state.json', 'final_response.md', 'coordinator_plan.md', 'neoag-sliding-run/sliding_run_summary.json', 'neoag-sliding-run/run-full/scoring/ranked_peptides.tsv', 'neoag-sliding-run/run-full/scoring/ranked_events.tsv', 'neoag-sliding-run/run-full/reports/evidence_report.html']
    files = []
    for rel in rels:
        p = outdir / rel
        if p.exists():
            files.append({'path': str(p), 'size': p.stat().st_size})
    return {'job': pub(job), 'outputs': files}

HTML = r"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>NeoAg Agent Internal Prototype</title>
  <style>
    body{font-family:Arial,sans-serif;margin:20px;background:#f6f7f9;color:#111827}
    section{background:white;padding:16px;margin:12px 0;border:1px solid #ddd;border-radius:8px}
    textarea,input,select{width:100%;box-sizing:border-box;margin:6px 0 12px;padding:8px;border:1px solid #cfd4dc;border-radius:4px}
    textarea{height:110px}
    button{padding:8px 12px;margin:4px;background:#155eef;color:white;border:0;border-radius:5px;cursor:pointer}
    button:disabled{background:#98a2b3;cursor:not-allowed}
    .danger{background:#b42318}
    .muted{color:#667085;font-size:13px}
    pre{background:#111;color:#eee;padding:12px;white-space:pre-wrap;max-height:520px;overflow:auto;border-radius:6px}
    .RUNNING{color:#155eef}.PASS{color:green}.FAIL,.STOPPED{color:#b42318}
    code{font-size:12px}
    .job-row{padding:6px 0;border-bottom:1px solid #eee}
  </style>
</head>
<body>
  <h2>NeoAg Agent <span class="muted">Internal prototype</span></h2>
  <section>
    <form id="jobForm">
      <p>自然语言任务</p>
      <textarea name="message">请描述要执行的任务，例如：检查输入文件并生成结果解读。</textarea>
      <p>上传 VCF / HLA 文件</p>
      <input id="vcfFile" name="vcf_upload" type="file" accept=".vcf,.gz,.txt,.tsv">
      <input id="hlaFile" name="hla_upload" type="file" accept=".txt,.tsv,.csv">
      <div class="muted">Internal prototype. Do not paste patient-identifying paths into shared demos. 可上传本地 VCF/HLA，也可填写服务器已有路径。</div>
      <p>服务器文件路径，每行一个</p>
      <textarea name="files" placeholder="/path/to/input.vcf.gz&#10;/path/to/hla.txt"></textarea>
      <p>输出目录</p>
      <input name="outdir" value="results/llm_agent_web">
      <p>Sample ID 过滤（可选，用于结果分析避免混入其他样本）</p>
      <input name="sample_id" placeholder="例如 SAMPLE001">
      <p>Provider / Model / Mode / API key env</p>
      <select name="llm_provider"><option>rule</option><option>deepseek</option><option>openai-compatible</option><option>vllm</option></select>
      <input name="model" placeholder="optional model name; leave blank for rule provider">
      <select name="mode"><option>execute-with-approval</option><option>execute-safe</option><option>plan</option></select>
      <input name="api_key_env" placeholder="optional API key env; leave blank for rule provider">
      <button id="runBtn" type="button">运行</button>
      <span id="notice" class="muted"></span>
    </form>
  </section>
  <section>
    <button id="refreshBtn" type="button">刷新任务</button>
    <button id="stopBtn" class="danger" type="button">停止当前</button>
    <h3>任务</h3>
    <div id="jobs">加载中...</div>
  </section>
  <section>
    <h3>最终回答</h3>
    <pre id="final">请选择或运行一个任务</pre>
    <h3>日志</h3>
    <div id="currentJob" class="muted">当前未选择任务</div>
    <pre id="log">请选择或运行一个任务</pre>
    <h3>输出</h3>
    <div id="outs" class="muted">暂无</div>
  </section>
<script>
(function(){
  var cur = null;
  var timer = null;
  function E(id){ return document.getElementById(id); }
  function setNotice(text){ E('notice').textContent = text || ''; }
  function scrollToLog(){ E('final').scrollIntoView({behavior:'smooth', block:'start'}); }
  async function api(path, opts){
    var r = await fetch(path, opts || {});
    var ct = r.headers.get('content-type') || '';
    var body = ct.indexOf('json') >= 0 ? await r.json() : await r.text();
    if(!r.ok){ throw new Error(typeof body === 'string' ? body : JSON.stringify(body)); }
    return body;
  }
  function lines(s){ return String(s || '').split("\n").map(function(x){ return x.trim(); }).filter(Boolean); }
  function makeUploadId(){ return new Date().toISOString().replace(/[^0-9]/g,'').slice(0,14) + '_' + Math.random().toString(16).slice(2,10); }
  async function uploadOne(file, uploadId){
    if(!file){ return null; }
    var url = '/api/upload?upload_id=' + encodeURIComponent(uploadId) + '&filename=' + encodeURIComponent(file.name);
    var r = await fetch(url, {method:'POST', headers:{'Content-Type':'application/octet-stream'}, body:file});
    var body = await r.json().catch(function(){ return {}; });
    if(!r.ok){ throw new Error(body.detail || JSON.stringify(body)); }
    return body.path;
  }
  async function refresh(){
    var d = await api('/api/jobs');
    E('jobs').innerHTML = d.jobs.map(function(j){
      return '<div class="job-row"><b class="'+j.status+'">'+j.status+'</b> <code>'+j.job_id+'</code> '+j.outdir+' <button type="button" data-job="'+j.job_id+'">查看</button></div>';
    }).join('') || '暂无任务';
    Array.prototype.forEach.call(E('jobs').querySelectorAll('button[data-job]'), function(b){
      b.addEventListener('click', function(){ selectJob(b.getAttribute('data-job')); });
    });
    if(!cur && d.jobs.length){
      var firstRunning = d.jobs.find(function(j){ return j.status === 'RUNNING'; }) || d.jobs[0];
      cur = firstRunning.job_id;
      E('currentJob').textContent = '当前任务: ' + cur + ' / ' + firstRunning.status;
    }
  }
  async function selectJob(id){
    cur = id;
    setNotice(' 正在查看 ' + id);
    E('currentJob').textContent = '当前任务: ' + id;
    if(timer){ clearInterval(timer); }
    timer = setInterval(showCurrent, 5000);
    await showCurrent();
    scrollToLog();
  }
  async function showCurrent(){
    if(!cur){ return; }
    var job = await api('/api/jobs/' + encodeURIComponent(cur));
    E('currentJob').textContent = '当前任务: ' + cur + ' / ' + job.status + ' / PID ' + (job.pid || '-');
    var finalText = await api('/api/jobs/' + encodeURIComponent(cur) + '/final');
    E('final').textContent = finalText && finalText.trim() ? finalText : '最终回答暂未生成。';
    var text = await api('/api/jobs/' + encodeURIComponent(cur) + '/log?tail=500');
    if(!text || !text.trim()){
      text = '任务已选择，但日志暂时为空。\n状态: ' + job.status + '\n日志文件: ' + (job.log_path || '-') + '\n输出目录: ' + (job.outdir || '-') + '\n命令: ' + (job.command || []).join(' ');
    }
    E('log').textContent = text;
    var o = await api('/api/jobs/' + encodeURIComponent(cur) + '/outputs');
    E('outs').innerHTML = o.outputs.map(function(x){ return '<div><code>'+x.path+'</code> '+x.size+' bytes</div>'; }).join('') || '暂无';
    await refresh();
  }
  async function stopCurrent(){
    if(!cur){ setNotice(' 请先选择任务'); return; }
    await api('/api/jobs/' + encodeURIComponent(cur) + '/stop', {method:'POST'});
    setNotice(' 已发送停止请求');
    await showCurrent();
  }
  async function runJob(){
    var btn = E('runBtn');
    btn.disabled = true;
    setNotice(' 正在上传/提交...');
    try{
      var fd = new FormData(E('jobForm'));
      var uploadId = makeUploadId();
      var files = lines(fd.get('files'));
      var vcfPath = await uploadOne(E('vcfFile').files[0], uploadId);
      var hlaPath = await uploadOne(E('hlaFile').files[0], uploadId);
      if(vcfPath){ files.push(vcfPath); }
      if(hlaPath){ files.push(hlaPath); }
      var body = {
        message: fd.get('message'),
        files: files,
        outdir: fd.get('outdir'),
        llm_provider: fd.get('llm_provider'),
        model: fd.get('model'),
        mode: fd.get('mode'),
        api_key_env: fd.get('api_key_env'),
        sample_id: fd.get('sample_id') || null,
        allow_high_risk: true
      };
      var j = await api('/api/jobs', {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify(body)});
      cur = j.job_id;
      E('currentJob').textContent = '当前任务: ' + j.job_id + ' / RUNNING';
      E('final').textContent = '新任务已提交，最终回答生成中。';
      E('log').textContent = '任务已提交，正在等待日志输出。';
      E('outs').textContent = '暂无';
      setNotice(' 已提交 ' + j.job_id);
      await refresh();
      await selectJob(j.job_id);
    }catch(err){
      setNotice(' 提交失败: ' + err.message);
    }finally{
      btn.disabled = false;
    }
  }
  E('jobForm').addEventListener('submit', function(e){ e.preventDefault(); runJob(); });
  E('runBtn').addEventListener('click', runJob);
  E('refreshBtn').addEventListener('click', refresh);
  E('stopBtn').addEventListener('click', stopCurrent);
  window.refresh = refresh;
  window.selectJob = selectJob;
  refresh().catch(function(err){ E('jobs').textContent = '加载失败: ' + err.message; });
})();
</script>
</body>
</html>"""
