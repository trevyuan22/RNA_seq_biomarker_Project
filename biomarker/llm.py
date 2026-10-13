"""Optional AI interpretation. Runs only if `anthropic` is installed AND ANTHROPIC_API_KEY is set.
Only a compact summary (gene names, pathway names, metrics) is sent - never the raw data."""
import json
import os

SYSTEM = ("You are a careful bioinformatics assistant. You will receive a JSON summary of an RNA-seq "
          "biomarker analysis. Write a short plain-language interpretation (max ~350 words) for a "
          "non-specialist: what each disease signature seems to involve, which pathways stand out, and "
          "how trustworthy the classifier looks. Use ONLY the genes, pathways and numbers provided; do not "
          "invent findings or cite papers. Clearly flag limitations. End with one sentence stating this "
          "is research output, not a medical diagnosis.")


def interpret(summary):
    if not os.environ.get("ANTHROPIC_API_KEY"):
        return None, "AI interpretation skipped: ANTHROPIC_API_KEY is not set."
    try:
        import anthropic
    except ImportError:
        return None, "AI interpretation skipped: `pip install anthropic` to enable it."
    try:
        client = anthropic.Anthropic()
        msg = client.messages.create(
            model=os.environ.get("CLAUDE_MODEL", "claude-sonnet-5-5"), max_tokens=1200, system=SYSTEM,
            messages=[{"role": "user", "content": json.dumps(summary, default=str)[:20000]}])
        return "".join(b.text for b in msg.content if getattr(b, "text", None)), ""
    except Exception as e:
        return None, f"AI interpretation failed: {type(e).__name__}"
